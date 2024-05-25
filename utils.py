import torch
import numpy as np
devict = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
def evaluate_model(model, test_loader, criterion):
    """
    Evaluate the performance of a given model.
    Args:
    model -- the PyTorch model to evaluate
    test_loader -- DataLoader for the test data
    criterion -- the loss function used
    Returns:
    Average test loss
    """
    model.eval()  # Set model to evaluation mode
    test_loss = 0.0
    total_mae_comp = 0.0
    total_mae_cont = 0.0
    total_count = 0.0
    with torch.no_grad():
        for data in test_loader:
            inputs, labels_comp, labels_cont = data
            inputs = inputs.to(device)
            labels_comp = labels_comp.to(device).unsqueeze(1)
            labels_cont = labels_cont.to(device).unsqueeze(1)
            outputs_comp, outputs_cont = model(inputs)
            loss_comp = criterion(outputs_comp, labels_comp)
            loss_cont = criterion(outputs_cont, labels_cont)
            test_loss += loss_comp.item() + loss_cont.item()
            mae_comp = torch.abs(outputs_comp - labels_comp).sum()  # Calculate MAE
            total_mae_comp += mae_comp.item()
            mae_cont = torch.abs(outputs_cont - labels_cont).sum()  # Calculate MAE
            total_mae_cont += mae_cont.item()
            total_count += labels_comp.size(0)  # Update count
    avg_mae_comp = total_mae_comp / total_count  # Calculate average MAE
    print('Test MAE_comp: %.3f' % avg_mae_comp)
    avg_mae_cont = total_mae_cont / total_count  # Calculate average MAE
    print('Test MAE_cont: %.3f' % avg_mae_cont)
    avg_test_loss = test_loss / len(test_loader)
    print('Test loss: %.3f' % avg_test_loss)
    return avg_test_loss, avg_mae_comp, avg_mae_cont