#train.py

import numpy as np
from scipy.sparse import csr_matrix
import scipy
import torch
from model import ResidualBlock, ResNet
from torch.utils.data import Dataset, DataLoader
import numpy as np
import torch
import time
import torch.nn as nn
import torchvision.transforms as transforms
# 检查CUDA是否可用，并设置device变量
device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
print(f"Using device: {device}")

#读取npy数据
train_x = scipy.sparse.load_npz('features/train_new_lineages_comp_scaled_vectors.npz').toarray()
train_y_comp= np.load('features/train_new_lineages_comp_scaled_labels_comp.npy',allow_pickle=True)/100
train_y_cont= np.load('features/train_new_lineages_comp_scaled_labels_cont.npy',allow_pickle=True)/100

#加载测试数据集
test_x = scipy.sparse.load_npz('features/test_pate_scaled_vectors.npz').toarray()
test_y_comp= np.load('features/test_pate_scaled_labels_comp.npy',allow_pickle=True)/100
test_y_cont= np.load('features/test_pate_scaled_labels_cont.npy',allow_pickle=True)/100
#将数据转换为张量
zero_padded_test_x=np.zeros((5996, 20164))
zero_padded_test_x[:, :20021] = test_x

target_shape = (1657133, 20164)
#创建一个新的形状为 (1337858, 20164) 的零数组
zero_padded_train_x = np.zeros(target_shape)
# 将 train_x 的值复制到新数组中
zero_padded_train_x[:, :20021] = train_x
del train_x, test_x
# augmentations = transforms.Compose([
#     transforms.RandomCrop(128),
#     transforms.RandomHorizontalFlip(),
#     transforms.RandomRotation(10),
#     #transforms.ColorJitter(brightness=0.1, contrast=0.1, saturation=0.1, hue=0.1),
#     transforms.RandomVerticalFlip(p=0.5),
#     transforms.RandomAffine(degrees=15, translate=(0.1, 0.1), scale=(0.9, 1.1), shear=10),
#     transforms.Normalize(mean=[0.5, 0.5, 0.5], std=[0.5, 0.5, 0.5])
#     # ... 其他需要的转换
# ])
class AddGaussianNoise(object):
    def __init__(self, mean=0., std=1.):
        self.mean = mean
        self.std = std

    def __call__(self, tensor):
        noise = torch.randn(tensor.size()) * self.std + self.mean
        return tensor + noise
class MetagenomeDataset(Dataset):
    def __init__(self, zero_padded_train_x, train_y_comp, train_y_cont, noise_transform=None):
        self.data = zero_padded_train_x
        self.label_comp = train_y_comp
        self.label2_cont = train_y_cont
        self.noise_transform = noise_transform

    def __len__(self):
        return len(self.data)

    def __getitem__(self, idx):
        image_shape = (142, 142)
        sample = self.data[idx]
        sample = sample.reshape(image_shape)

        # 添加噪音
        if self.noise_transform:
            sample = self.noise_transform(torch.Tensor(sample))

        label_comp = self.label_comp[idx]
        label_cont = self.label2_cont[idx]

        # 将sample、label_comp和label_cont转换为PyTorch张量
        sample = (torch.Tensor(sample)).unsqueeze(0)
        label_comp = torch.Tensor(np.array(label_comp))
        label_cont = torch.Tensor(np.array(label_cont))
        return sample, label_comp, label_cont
    
def evaluate_model(model, test_loader, criterion):
    """
    评估给定模型的性能。
    参数:
    model -- 要评估的PyTorch模型
    test_loader -- 测试数据的DataLoader
    criterion -- 使用的损失函数
    返回:
    测试集上的平均损失
    """
    model.eval()  # 设置模型为评估模式
    test_loss = 0.0
    total_mae_comp = 0.0
    total_mae_cont=0.0
    total_count = 0.0
    with torch.no_grad():
        for data in test_loader:
            inputs, labels_comp, labels_cont = data
            inputs = inputs.to(device)
            labels_comp = labels_comp.to(device).unsqueeze(1)
            labels_cont = labels_cont.to(device).unsqueeze(1) 
            outputs_comp, outputs_cont = model(inputs)
            loss_comp = criterion(outputs_comp, labels_comp)  # 假设您想预测的是 labels_comp
            loss_cont = criterion(outputs_cont, labels_cont)  # 假设您想预测的是 labels_cont
            test_loss += loss_comp.item()+loss_cont.item()
            mae_comp = torch.abs(outputs_comp - labels_comp).sum()  # 计算MAE
            total_mae_comp += mae_comp.item()
            mae_cont = torch.abs(outputs_cont - labels_cont).sum()  # 计算MAE
            total_mae_cont += mae_cont.item()
            total_count += labels_comp.size(0)  # 更新计数
    avg_mae_comp = total_mae_comp / total_count  # 计算平均MAE
    print('Test MAE_comp: %.3f' % avg_mae_comp)
    avg_mae_cont = total_mae_cont / total_count  # 计算平均MAE
    print('Test MAE_cont: %.3f' % avg_mae_cont)
    avg_test_loss = test_loss / len(test_loader)
    print('Test loss: %.3f' % avg_test_loss)
    return avg_test_loss, avg_mae_comp, avg_mae_cont
noise_transform = AddGaussianNoise(0, 0.2)
Metagenome_dataset= MetagenomeDataset(zero_padded_train_x, train_y_comp, train_y_cont)
train_loader = DataLoader(Metagenome_dataset, batch_size=512, shuffle=True)
Metagenome_test=MetagenomeDataset(zero_padded_test_x,test_y_comp,test_y_cont)
test_loader = DataLoader(Metagenome_test, batch_size=512, shuffle=False)
model = ResNet(ResidualBlock, [2, 2, 2, 2]).to(device)
import torch.optim as optim
criterion = nn.MSELoss()
#criterion = nn.SmoothL1Loss()
##########################################################################
optimizer = optim.Adam(model.parameters(), lr=0.0001, weight_decay=1e-5)
num_epochs = 10000  # 或您选择的任何epoch数
best_avg_mae=0.12
for epoch in range(num_epochs):
    model.train()  # 设置模型为训练模式
    running_loss = 0.0
    start_time = time.time()
    for i, data in enumerate(train_loader):
        # 获取输入数据
        inputs, labels_comp, labels_cont = data
        inputs = inputs.to(device)
        labels_comp = labels_comp.to(device).unsqueeze(1) 
        labels_cont = labels_cont.to(device).unsqueeze(1) 
        # 梯度清零
        optimizer.zero_grad()
        # 前向传播
        outputs_comp, outputs_cont= model(inputs)
        # 计算损失
        comp_loss = criterion(outputs_comp, labels_comp)  # 假设您想预测的是 labels_comp
        cont_loss = criterion(outputs_cont, labels_cont)  # 假设您想预测的是 labels_cont
        loss=0.5*comp_loss+cont_loss
        # 反向传播
        loss.backward()
        # 优化
        optimizer.step()
        running_loss += loss.item()
        if i % 1000 == 999:    # 每1000个mini-batches打印一次
            print('[%d, %5d] loss: %.3f' %
                  (epoch + 1, i + 1, running_loss / 1000))
            running_loss = 0.0
    end_time = time.time()
    training_time_seconds = (end_time - start_time)/3600
    # 打印训练时间
    print(f"训练时间: {training_time_seconds:.2f} 小时")
    avg_test_loss, avg_mae_comp, avg_mae_cont=evaluate_model(model, test_loader, criterion)
    if avg_mae_cont<best_avg_mae:
        best_avg_mae=avg_mae_cont
        model_save_name = "models/deepcheck_pate/deepcheck_"+str(avg_mae_cont)+"new.pt"
        torch.save(model.state_dict(), model_save_name)
        print("模型已保存")