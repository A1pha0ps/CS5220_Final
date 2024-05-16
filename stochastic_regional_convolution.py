import torch
import torch.nn as nn
import torch.nn.functional as F
import numpy as np

class StochasticRegionalConvolution(nn.Module):
    """
    Custom PyTorch layer for Stochastic Regional Convolutions (SRCs).

    Args:
        in_channels (int): Number of input channels.
        out_channels (int): Number of output channels.
        kernel_size (int): Size of the convolution kernel.
        region_size (int): Size of the region to be convolved.
        T (int): Number of regions to randomly select and convolve.

    Forward pass:
        - Selects T random regions from the input.
        - Applies convolution to each selected region.
        - Aggregates the results using a weighted sum.
    """
    def __init__(self, in_channels, out_channels, kernel_size, region_size, T):
        super(StochasticRegionalConvolution, self).__init__()
        self.in_channels = in_channels
        self.out_channels = out_channels
        self.kernel_size = kernel_size
        self.region_size = region_size
        self.T = T
        
        # Define a convolutional layer
        self.conv = nn.Conv2d(in_channels, out_channels, kernel_size, bias=False)

    def forward(self, x):
        B, C, H, W = x.size()
        k = self.kernel_size
        r = self.region_size
        n_r = (H // r) * (W // r)
        
        # Initialize output tensor
        output = torch.zeros_like(x)
        
        # Perform convolution on T random regions
        for _ in range(self.T):
            h_i = np.random.randint(0, H - r + 1)
            w_i = np.random.randint(0, W - r + 1)
            x_i = x[:, :, h_i:h_i + r, w_i:w_i + r]
            c_i = self.conv(x_i)
            lambda_i = np.random.rand()
            output[:, :, h_i:h_i + r - k + 1, w_i:w_i + r - k + 1] += lambda_i * c_i
        
        # Normalize the output
        output /= n_r
        return output

class SRC_CNN(nn.Module):
    """
    A simple CNN model that integrates the custom Stochastic Regional Convolution layer.

    Architecture:
        - StochasticRegionalConvolution layer
        - MaxPool layer
        - Conv2d layer
        - Fully connected layers
    """
    def __init__(self):
        super(SRC_CNN, self).__init__()
        self.src1 = StochasticRegionalConvolution(3, 16, 3, 10, 50)
        self.pool = nn.MaxPool2d(2, 2)
        self.conv2 = nn.Conv2d(16, 32, 3)
        self.fc1 = nn.Linear(32 * 6 * 6, 120)
        self.fc2 = nn.Linear(120, 84)
        self.fc3 = nn.Linear(84, 10)

    def forward(self, x):
        x = self.pool(F.relu(self.src1(x)))
        x = self.pool(F.relu(self.conv2(x)))
        x = x.view(-1, 32 * 6 * 6)
        x = F.relu(self.fc1(x))
        x = F.relu(self.fc2(x))
        x = self.fc3(x)
        return x


