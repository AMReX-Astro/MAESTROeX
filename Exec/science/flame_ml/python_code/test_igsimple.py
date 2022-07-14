# Python script to train an ML model for ignition_simple network
import os
import shutil
import random
import torch
import torch.nn as nn
import torch.optim as optim

from train import NuclearReactionML
from networks import *
from losses import logX_loss_noenuc, loss_mass_fraction_log_conserv, loss_wexp_noenuc, loss_mass_fraction_conserv


## FILE SETUP

# this is the path to your data files
data_path = '../data/'

# These is the input/output prefix of your datafile names.
input_prefix = 'react_inputs_*'
output_prefix = 'react_outputs_*'

# Plotfile prefixes, used for visualization purposes.
plotfile_prefix = 'flame_*'

# By default, this package will save your model, logs of the training and testing data during training,
# and plots to a directory. Here you specify that directory.
output_dir = 'testing123/'

# Check to see if the output directory already exists
if os.path.exists(output_dir):
    hash = random.getrandbits(32)
    new_output_dir = output_dir[:-1] + f"_{hex(hash)}/"
    print(f"Replacing output directory {output_dir} with {new_output_dir}")
    #os.rename(output_dir, new_output_dir)
    output_dir = new_output_dir

# The log file. Everything that is printed during training also goes into this file in case something
# gets interrupted.
log_file = output_dir + "log.txt"


## MODEL SETUP

nrml = NuclearReactionML(data_path, input_prefix, output_prefix, plotfile_prefix,
                         output_dir, log_file, DEBUG_MODE=False, DO_PLOTTING=True,
                         SAVE_MODEL=True, DO_HYPER_OPTIMIZATION=False, LOG_MODE=True)

device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
print(device)
if torch.cuda.is_available():
    print(f"Using {torch.cuda.device_count()} GPUs!")


def criterion(pred, target): 
    loss1 = 10*loss_mass_fraction_conserv(pred, target, nnuc=2)
    #loss1 = 10*loss_mass_fraction_log_conserv(1.0, pred, target, nnuc=2)
    loss2 = loss_wexp_noenuc(pred, target, nnuc=2, offset=10.0)
    #loss2 = logX_loss_noenuc(pred, target, nnuc=2)
    
    L = nn.MSELoss()
    F = nn.L1Loss()
    enuc_fac = 1
    loss_enuc = L(pred[:,2], target[:,2]) + enuc_fac * F(torch.sign(pred[:,2]), torch.sign(target[:,2]))
    return loss_enuc + loss2 + loss1

    
## TRAIN MODEL

model_spec = ResNet(4, 16, 32, 32, 16, 2, relu=True)
model_spec.to(device)
    
model_enuc = Combine_Net3(4, 32, 16, 8, 8, 8, 8, 8, 8, 16, 32, 1, relu=False)
model_enuc.to(device)

model = Combine2Models(model_spec, model_enuc)
print(model)
# get model to cuda if possible
model.to(device)

num_epochs = 10

optimizer = optim.Adam(model.parameters(), lr=1e-5)

nrml.train(model, optimizer, num_epochs, criterion)

# Plot at the end
nrml.plot()


## SAVE MODEL

device_cpu = torch.device('cpu')

# reload model onto cpu
model_spec_cpu = ResNet(4, 16, 32, 32, 16, 2, relu=True)
model_spec_cpu.to(device_cpu)
model_enuc_cpu = Combine_Net3(4, 32, 16, 8, 8, 8, 8, 8, 8, 16, 32, 1, relu=False)
model_enuc_cpu.to(device_cpu)
model_cpu = Combine2Models(model_spec_cpu, model_enuc_cpu)
print("Loading model onto CPU...")

model_cpu.load_state_dict(torch.load(output_dir + "my_model.pt", map_location=device_cpu))
print(model_cpu)

# convert to torch script
print("Saving model...")
net_module = torch.jit.script(model_cpu)
net_module.save(output_dir + "ts_model.pt")
print(net_module.code)
