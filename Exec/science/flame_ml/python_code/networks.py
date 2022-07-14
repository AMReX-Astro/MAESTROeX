import torch.nn as nn
import torch

class Combine2Models(nn.Module):
    def __init__(self, ModelSpec, ModelEnuc):
        super().__init__()
        self.ModelSpec = ModelSpec
        self.ModelEnuc = ModelEnuc

    def forward(self, x):
        xspec = self.ModelSpec(x)
        xenuc = self.ModelEnuc(x)
        return torch.cat((xspec, xenuc), dim=1)


# Nets inspired by ResNet
class ResNet(nn.Module):
    def __init__(self, input_size, h1, h2, h3, h4, output_size,
                 relu=False, tanh=True):
        super().__init__()
        self.relu = relu
        self.fc1 = nn.Linear(input_size, h1)
        self.ac1 = nn.Tanh() if tanh else nn.CELU()
        self.fc2 = nn.Linear(h1, h2)
        self.ac2 = nn.Tanh() if tanh else nn.CELU()
        self.fc3 = nn.Linear(h2, h3)
        self.ac3 = nn.Tanh() if tanh else nn.CELU()
        self.fc4 = nn.Linear(h3, h4)
        self.ac4 = nn.Tanh() if tanh else nn.CELU()
        self.fc5 = nn.Linear(h4, output_size)
        self.ac5 = nn.ReLU()
        
        # layers between non-consecutive layers 
        self.fc0to3 = nn.Linear(input_size, h3)
        self.fc2to5 = nn.Linear(h2, output_size)

        if self.relu:
            # initialize final layer with He weights that work better with ReLU activation
            nn.init.kaiming_normal_(self.fc5.weight, nonlinearity='relu')
            nn.init.kaiming_normal_(self.fc2to5.weight, nonlinearity='relu')

    def forward(self, x):
        x1 = self.ac1(self.fc1(x))
        x2 = self.ac2(self.fc2(x1))
        x3 = self.ac3(self.fc3(x2) + self.fc0to3(x))
        x4 = self.ac4(self.fc4(x3))
        x5 = self.fc5(x4) + self.fc2to5(x2)
        if self.relu:
            x5 = self.ac5(x5)
        return x5

class Combine_Net3(nn.Module):
    def __init__(self, input_size, h1, h2, h3, h4, h5, h6, h7, h8, h9, h10, output_size,
                 relu=False, tanh=True):
        super().__init__()
        self.relu = relu
        self.fc1 = nn.Linear(input_size, h1)
        self.ac1 = nn.Tanh() if tanh else nn.CELU()
        self.fc2 = nn.Linear(h1, h2)
        self.ac2 = nn.Tanh() if tanh else nn.CELU()
        self.fc3 = nn.Linear(h2, h3)
        self.ac3 = nn.Tanh() if tanh else nn.CELU()
        self.fc4 = nn.Linear(h3, h4)
        self.ac4 = nn.Tanh() if tanh else nn.CELU()
        self.fc5 = nn.Linear(h4, h5)
        self.ac5 = nn.Tanh() if tanh else nn.CELU()
        self.fc6 = nn.Linear(h5, h6)
        self.ac6 = nn.Tanh() if tanh else nn.CELU()
        self.fc7 = nn.Linear(h6, h7)
        self.ac7 = nn.Tanh() if tanh else nn.CELU()
        self.fc8 = nn.Linear(h7, h8)
        self.ac8 = nn.Tanh() if tanh else nn.CELU()
        self.fc9 = nn.Linear(h8, h9)
        self.ac9 = nn.Tanh() if tanh else nn.CELU()
        self.fc10 = nn.Linear(h9, h10)
        self.ac10 = nn.Tanh() if tanh else nn.CELU()
        self.fc11 = nn.Linear(h10, output_size)
        self.ac11 = nn.ReLU()  # same as removing negative values
        
        # Resnet connections 
        self.fc1to4 = nn.Linear(h1, h4)
        self.fc3to6 = nn.Linear(h3, h6)
        self.fc5to8 = nn.Linear(h5, h8)
        self.fc7to10 = nn.Linear(h7, h10)
        
        # U-net-like connections
        self.fc1to10 = nn.Linear(h1, h10)
        self.fcio = nn.Linear(input_size, output_size)

        if self.relu:
            # initialize final layer with He weights that work better with ReLU activation
            nn.init.kaiming_normal_(self.fc11.weight, nonlinearity='relu')
            nn.init.kaiming_normal_(self.fcio.weight, nonlinearity='relu')

    def forward(self, x):
        x1 = self.ac1(self.fc1(x))
        x2 = self.ac2(self.fc2(x1))
        x3 = self.ac3(self.fc3(x2))
        x4 = self.ac4(self.fc4(x3) + self.fc1to4(x1))
        x5 = self.ac5(self.fc5(x4))
        x6 = self.ac6(self.fc6(x5) + self.fc3to6(x3))
        x7 = self.ac7(self.fc7(x6))
        x8 = self.ac8(self.fc8(x7) + self.fc5to8(x5))
        x9 = self.ac9(self.fc9(x8))
        x10 = self.ac10(self.fc10(x9) + self.fc7to10(x7) + self.fc1to10(x1))
        x11 = self.fc11(x10) + self.fcio(x)
        if self.relu:
            x11 = self.ac11(x11)
        return x11
