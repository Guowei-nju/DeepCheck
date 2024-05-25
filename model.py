import torch
import torch.nn as nn
import torch.nn.functional as F

# 定义一个残差块
class ResidualBlock(nn.Module):
    def __init__(self, in_channels, out_channels, stride=1, downsample=None):
        super(ResidualBlock, self).__init__()
        # 第一个卷积层
        self.conv1 = nn.Conv2d(in_channels, out_channels, kernel_size=3, stride=stride, padding=1, bias=False)
        self.bn1 = nn.BatchNorm2d(out_channels)
        self.relu = nn.ReLU(inplace=True)
        # 第二个卷积层
        self.conv2 = nn.Conv2d(out_channels, out_channels, kernel_size=3, padding=1, bias=False)
        self.bn2 = nn.BatchNorm2d(out_channels)
        self.downsample = downsample

    def forward(self, x):
        identity = x
        out = self.conv1(x)
        out = self.bn1(out)
        out = self.relu(out)
        out = self.conv2(out)
        out = self.bn2(out)

        if self.downsample is not None:
            identity = self.downsample(x)

        out += identity
        out = self.relu(out)
        return out
class SelfAttention(nn.Module):
    def __init__(self, in_channels):
        super(SelfAttention, self).__init__()
        # 定义自注意力的键、查询和值
        self.query_conv = nn.Conv2d(in_channels, in_channels // 8, kernel_size=1)
        self.key_conv = nn.Conv2d(in_channels, in_channels // 8, kernel_size=1)
        self.value_conv = nn.Conv2d(in_channels, in_channels, kernel_size=1)
        self.softmax = nn.Softmax(dim=-1)
        
    def forward(self, x):
        # [Batch, Channel, Height, Width]
        batch_size, C, height, width = x.size()
        # 查询
        query = self.query_conv(x).view(batch_size, -1, height * width).permute(0, 2, 1)
        # 键
        key = self.key_conv(x).view(batch_size, -1, height * width)
        # 注意力分数
        attention = self.softmax(torch.bmm(query, key))
        # 值
        value = self.value_conv(x).view(batch_size, -1, height * width)
        # 输出
        out = torch.bmm(value, attention.permute(0, 2, 1))
        out = out.view(batch_size, C, height, width)
        
        return out
# 定义一个ResNet模型
class ResNet(nn.Module):
    def __init__(self, block, layers, num_classes=1):
        super(ResNet, self).__init__()
        self.in_channels = 64
        # 第一个卷积层
        self.conv1 = nn.Conv2d(1, 64, kernel_size=7, stride=2, padding=3, bias=False)
        self.bn1 = nn.BatchNorm2d(64)
        self.relu = nn.ReLU(inplace=True)
        self.maxpool = nn.MaxPool2d(kernel_size=3, stride=2, padding=1)
        # 四个残差层
        self.layer1 = self._make_layer(block, 64, layers[0])
        self.layer2 = self._make_layer(block, 128, layers[1], stride=2)
        self.layer3 = self._make_layer(block, 256, layers[2], stride=2)
        self.layer4 = self._make_layer(block, 512, layers[3], stride=2)
        self.attention = SelfAttention(512)
        self.avgpool = nn.AdaptiveAvgPool2d((1, 1))
        # 全连接层
        self.dropout = nn.Dropout(0.5)
        self.fc = nn.Linear(512, 100)
        self.fc1 = nn.Linear(100, num_classes)
        self.fc2 = nn.Linear(100, num_classes)


    def _make_layer(self, block, out_channels, blocks, stride=1):
        downsample = None
        if stride != 1 or self.in_channels != out_channels:
            downsample = nn.Sequential(
                nn.Conv2d(self.in_channels, out_channels, kernel_size=1, stride=stride, bias=False),
                nn.BatchNorm2d(out_channels))
        layers = []
        layers.append(block(self.in_channels, out_channels, stride, downsample))
        self.in_channels = out_channels
        for _ in range(1, blocks):
            layers.append(block(out_channels, out_channels))
        return nn.Sequential(*layers)

    def forward(self, x):
        x = self.conv1(x)
        x = self.bn1(x)
        x = self.relu(x)
        x = self.maxpool(x)
        x = self.layer1(x)
        x = self.layer2(x)
        x = self.layer3(x)
        x = self.layer4(x)

        # 应用自注意力模块
        x = self.attention(x)
        x = self.avgpool(x)
        x = torch.flatten(x, 1)
        x=self.dropout(x)
        x=self.fc(x)
        x_comp = self.fc1(x)
        x_cont = self.fc2(x)
        return x_comp
#        x = self.fc(x)
#        return x
class ResNet_cover_conv(nn.Module):
    def __init__(self, block, layers, num_classes=1):
        super(ResNet_cover_conv, self).__init__()
        self.in_channels = 64
        # 第一个卷积层
        self.conv1 = nn.Conv2d(1, 64, kernel_size=7, stride=2, padding=3, bias=False)
        self.bn1 = nn.BatchNorm2d(64)
        self.relu = nn.ReLU(inplace=True)
        self.maxpool = nn.MaxPool2d(kernel_size=3, stride=2, padding=1)
        # 四个残差层
        self.layer1 = self._make_layer(block, 64, layers[0])
        self.layer2 = self._make_layer(block, 128, layers[1], stride=2)
        self.layer3 = self._make_layer(block, 256, layers[2], stride=2)
        self.layer4 = self._make_layer(block, 512, layers[3], stride=2)
        self.layer5 = self._make_layer(block, 512, layers[4], stride=2)
        self.attention1 = SelfAttention(512)
        self.attention2 = SelfAttention(512)
        self.avgpool = nn.AdaptiveAvgPool2d((1, 1))
        # 全连接层
        self.dropout = nn.Dropout(0.5)
        self.fc_comp = nn.Linear(512, 100)
        self.fc_cont = nn.Linear(512, 100)
        self.fc = nn.Linear(512, 100)
        self.fc1 = nn.Linear(100, num_classes)
        self.fc2 = nn.Linear(100, num_classes)


    def _make_layer(self, block, out_channels, blocks, stride=1):
        downsample = None
        if stride != 1 or self.in_channels != out_channels:
            downsample = nn.Sequential(
                nn.Conv2d(self.in_channels, out_channels, kernel_size=1, stride=stride, bias=False),
                nn.BatchNorm2d(out_channels))
        layers = []
        layers.append(block(self.in_channels, out_channels, stride, downsample))
        self.in_channels = out_channels
        for _ in range(1, blocks):
            layers.append(block(out_channels, out_channels))
        return nn.Sequential(*layers)

    def forward(self, x):
        x = self.conv1(x)
        x = self.bn1(x)
        x = self.relu(x)
        x = self.maxpool(x)
        x = self.layer1(x)
        x = self.layer2(x)
        x = self.layer3(x)
        x1 = self.layer4(x)
        x2 = self.layer5(x)
        # 应用自注意力模块
        x1 = self.attention1(x1)
        x1 = self.avgpool(x1)
        x1 = torch.flatten(x1, 1)
        x1=self.dropout(x1)
        x1=self.fc_comp(x1)
        x_comp = self.fc1(x1)
    ################
        x2 = self.attention2(x2)
        x2 = self.avgpool(x2)
        x2 = torch.flatten(x2, 1)
        x2=self.dropout(x2)
        x2=self.fc_cont(x2)
        x_cont = self.fc2(x2)
        return x_comp,x_cont