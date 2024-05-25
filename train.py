import numpy as np
from scipy.sparse import load_npz
import torch
from model import ResidualBlock, ResNet
from torch.utils.data import Dataset, DataLoader
import torch.nn as nn
import torchvision.transforms as transforms
from sklearn.model_selection import train_test_split
import argparse
import time
import os
from utils import evaluate_model

# Set up argument parser
parser = argparse.ArgumentParser(description='Train a model on metagenomic data.')
parser.add_argument('--train_x', type=str, required=True, help='Path to the training feature file')
parser.add_argument('--train_y_comp', type=str, required=True, help='Path to the training comp labels file')
parser.add_argument('--train_y_cont', type=str, required=True, help='Path to the training cont labels file')
parser.add_argument('--test_x', type=str, required=False, help='Path to the testing feature file')
parser.add_argument('--test_y_comp', type=str, required=False, help='Path to the testing comp labels file')
parser.add_argument('--test_y_cont', type=str, required=False, help='Path to the testing cont labels file')
parser.add_argument('--batch_size', type=int, default=512, help='Batch size for training')
parser.add_argument('--num_epochs', type=int, default=10000, help='Number of epochs to train')
parser.add_argument('--learning_rate', type=float, default=0.0001, help='Learning rate')
parser.add_argument('--weight_decay', type=float, default=1e-5, help='Weight decay for optimizer')
args = parser.parse_args()

# Check if CUDA is available and set device variable
device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
print(f"Using device: {device}")

# Read npy data
train_x = load_npz(args.train_x).toarray()
train_y_comp = np.load(args.train_y_comp, allow_pickle=True) / 100
train_y_cont = np.load(args.train_y_cont, allow_pickle=True) / 100

# Check if test files are provided and load test data if available
if args.test_x and args.test_y_comp and args.test_y_cont and os.path.exists(args.test_x) and os.path.exists(args.test_y_comp) and os.path.exists(args.test_y_cont):
    test_x = load_npz(args.test_x).toarray()
    test_y_comp = np.load(args.test_y_comp, allow_pickle=True) / 100
    test_y_cont = np.load(args.test_y_cont, allow_pickle=True) / 100

    zero_padded_test_x = np.zeros((test_x.shape[0], 20164))
    zero_padded_test_x[:, :20021] = test_x
else:
    test_x, test_y_comp, test_y_cont = None, None, None

# Convert data to tensors
target_shape = (train_x.shape[0], 20164)
zero_padded_train_x = np.zeros(target_shape)
zero_padded_train_x[:, :20021] = train_x
del train_x

# Split the training data into training and validation sets
train_x, val_x, train_y_comp, val_y_comp, train_y_cont, val_y_cont = train_test_split(
    zero_padded_train_x, train_y_comp, train_y_cont, test_size=0.1, random_state=42
)

class AddGaussianNoise(object):
    def __init__(self, mean=0., std=1.):
        self.mean = mean
        self.std = std

    def __call__(self, tensor):
        noise = torch.randn(tensor.size()) * self.std + self.mean
        return tensor + noise

class MetagenomeDataset(Dataset):
    def __init__(self, data, label_comp, label_cont, noise_transform=None):
        self.data = data
        self.label_comp = label_comp
        self.label_cont = label_cont
        self.noise_transform = noise_transform

    def __len__(self):
        return len(self.data)

    def __getitem__(self, idx):
        image_shape = (142, 142)
        sample = self.data[idx]
        sample = sample.reshape(image_shape)

        # Add noise
        if self.noise_transform:
            sample = self.noise_transform(torch.Tensor(sample))

        label_comp = self.label_comp[idx]
        label_cont = self.label_cont[idx]

        # Convert sample, label_comp, and label_cont to PyTorch tensors
        sample = (torch.Tensor(sample)).unsqueeze(0)
        label_comp = torch.Tensor(np.array(label_comp))
        label_cont = torch.Tensor(np.array(label_cont))
        return sample, label_comp, label_cont

# Data loading
noise_transform = AddGaussianNoise(0, 0.2)
train_dataset = MetagenomeDataset(train_x, train_y_comp, train_y_cont)
val_dataset = MetagenomeDataset(val_x, val_y_comp, val_y_cont)

train_loader = DataLoader(train_dataset, batch_size=args.batch_size, shuffle=True)
val_loader = DataLoader(val_dataset, batch_size=args.batch_size, shuffle=False)

if test_x is not None:
    test_dataset = MetagenomeDataset(zero_padded_test_x, test_y_comp, test_y_cont)
    test_loader = DataLoader(test_dataset, batch_size=args.batch_size, shuffle=False)

# Model, loss function, and optimizer
model = ResNet(ResidualBlock, [2, 2, 2, 2]).to(device)
criterion = nn.MSELoss()
optimizer = torch.optim.Adam(model.parameters(), lr=args.learning_rate, weight_decay=args.weight_decay)
best_avg_mae=0
# Training loop
for epoch in range(args.num_epochs):
    model.train()  # Set model to training mode
    running_loss = 0.0
    start_time = time.time()
    for i, data in enumerate(train_loader):
        # Get input data
        inputs, labels_comp, labels_cont = data
        inputs = inputs.to(device)
        labels_comp = labels_comp.to(device).unsqueeze(1)
        labels_cont = labels_cont.to(device).unsqueeze(1)
        
        # Zero the parameter gradients
        optimizer.zero_grad()
        
        # Forward pass
        outputs_comp, outputs_cont = model(inputs)
        
        # Compute loss
        comp_loss = criterion(outputs_comp, labels_comp)
        cont_loss = criterion(outputs_cont, labels_cont)
        loss = 0.5 * comp_loss + cont_loss
        
        # Backward pass
        loss.backward()
        
        # Optimize
        optimizer.step()
        
        running_loss += loss.item()
        if i % 1000 == 999:  # Print every 1000 mini-batches
            print('[%d, %5d] loss: %.3f' % (epoch + 1, i + 1, running_loss / 1000))
            running_loss = 0.0
    end_time = time.time()
    training_time_seconds = (end_time - start_time) / 3600
    print(f"Training time: {training_time_seconds:.2f} hours")
    
    avg_val_loss, avg_val_mae_comp, avg_val_mae_cont = evaluate_model(model, val_loader, criterion)
    if avg_val_mae_cont < best_avg_mae:
        best_avg_mae = avg_val_mae_cont
        model_save_name = f"models/eepcheck_{avg_val_mae_cont:.4f}_new.pt"
        torch.save(model.state_dict(), model_save_name)
        print("Model saved")

# Finally evaluate the model on the test set if
if test_x is not None
    evaluate_model(model, test_loader, criterion)
