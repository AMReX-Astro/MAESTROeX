import yt
import numpy as np
import matplotlib.pyplot as plt
import matplotlib

import os
from glob import glob
import warnings
import sys
import pandas as pd
from datetime import datetime
import torch
import torch.nn as nn
import torch.nn.functional as F
import torch.optim as optim
import torch.utils.data
from torch.utils.data.dataset import Dataset
from torch.utils.data import DataLoader
from reactdataset_nodt import ReactDataset
from losses import component_loss_f, component_loss_f_L1
from plotting import plotting_standard


yt.funcs.mylog.setLevel(40) # Gets rid of all of the yt info text, only errors.
warnings.filterwarnings("ignore",category=matplotlib.cbook.mplDeprecation) #ignore plt depreciations
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")


class NuclearReactionML:


    def __init__(self, data_path, input_prefix, output_prefix, plotfile_prefix,
                 output_dir, log_file, DEBUG_MODE=True, DO_PLOTTING=True,
                 SAVE_MODEL=True, DO_HYPER_OPTIMIZATION=False, LOG_MODE=True):
                """
                data_path (string): this is the path to your data files

                input_prefix (string): This is the input prefix of your
                                       datafile names.

                output_prefix (string): This is the output prefix of your
                                        datafile names.

                plotfile_prefix (string): Plotfile prefixes, used for
                                          visualization purposes

                output_dir (string): By default, this package will save your model,
                                     logs of the training and testing data during
                                     training, and plots to a directory. Here
                                     you specify that directory. Must be a new
                                     directory or you will get an error to
                                     prevent the overwrite of data.

                log_file (string): The log file. Everything that is printed
                                   during training also goes into this file
                                   in case something gets interrupted.


                DEBUG_MODE (bool): This takes a small cut of the data (5 plotfiles)
                                   to train on. This is useful for debuggind since
                                   you won't have to deal with unwieldy amounts
                                   of data. Default=True

                DO_PLOTTING (bool): Whether to do the error plots or not.
                                    Saved in output_dir Default=True

                SAVE_MODEL (bool):  Whether to save the pytorch model file or
                                    not. Saved in output_dir Deafult=True

                DO_HYPER_OPTIMIZATION (bool): This still needs to be tweaked.
                """

                self.DO_PLOTTING = DO_PLOTTING
                self.output_dir = output_dir
                self.SAVE_MODEL = SAVE_MODEL
                self.DO_HYPER_OPTIMIZATION = DO_HYPER_OPTIMIZATION
                self.LOG_MODE = LOG_MODE

                if os.path.isdir(output_dir) and (len(os.listdir(output_dir)) != 0):
                    print(f"Directory {output_dir} exists and is not empty.")
                    print("Please change output_dir or remove the directory to prevent overwritting data.")
                    sys.exit()

                isdir = os.path.isdir(output_dir)
                if not isdir:
                    os.mkdir(output_dir)

                from tools import Logger
                self.logger = Logger(log_file)

                now = datetime.now()
                dt_string = now.strftime("%m/%d/%Y %H:%M:%S")
                self.logger.write(f"Model starting on : {dt_string}")
                self.logger.write(f"input_prefix {input_prefix}")
                self.logger.write(f"output_prefix {output_prefix}")
                self.logger.write(f"output_dir {output_dir}")
                self.logger.write(f"DEBUG_MODE {DEBUG_MODE}")
                self.logger.write(f"DO_PLOTTING {DO_PLOTTING}")
                #self.logger.write(f"log_loss_option {log_loss_option}")
                self.logger.write(f"DO_HYPER_OPTIMIZATION {DO_HYPER_OPTIMIZATION}")


                #LOADING DATA----------------------------------------------------------
                plotfiles = glob(data_path + plotfile_prefix)
                plotfiles = sorted(plotfiles)
                plotfiles = plotfiles[:-2] #cut after divuiter and initproj
                plotfiles = [plotfiles[-1]] + plotfiles[:-1] #move initdata to front.
                #make_movie(plotfiles, movie_name='enuc.mp4', var='enuc')

                react_data = ReactDataset(data_path, input_prefix, output_prefix, plotfile_prefix, DEBUG_MODE=DEBUG_MODE)
                self.nnuc = int(react_data.output_data.shape[1]/2 - 1)

                #Normalize density, temperature, and enuc
                dens_fac = torch.max(react_data.input_data[:, self.nnuc+1, :])
                temp_fac = torch.max(react_data.input_data[:, self.nnuc+2, :])
                enuc_fac = 1.1*torch.max(react_data.output_data[:, self.nnuc, :])
                react_data.input_data[:, self.nnuc+1, :]  = react_data.input_data[:, self.nnuc+1, :]/dens_fac
                react_data.input_data[:, self.nnuc+2, :]  = react_data.input_data[:, self.nnuc+2, :]/temp_fac
                react_data.output_data[:, self.nnuc, :] = react_data.output_data[:, self.nnuc, :]/enuc_fac
                
                #save these factors to a file
                arr = np.array([dens_fac.item(), temp_fac.item(), enuc_fac.item()])
                np.savetxt(self.output_dir + 'scaling_factors.txt', arr, header='Density, Temperature, Enuc factors (ordered)')

                if self.LOG_MODE:
                    #take 1/log of mass fractions of species 
                    react_data.input_data[:,1:self.nnuc+1,:] = -1.0/torch.log(react_data.input_data[:,1:self.nnuc+1,:])
                    react_data.output_data[:,:self.nnuc,:] = -1.0/torch.log(react_data.output_data[:,:self.nnuc,:])

                self.fields = [field for field in yt.load(react_data.output_files[0])._field_list]
                #truncate to mass fractions + enuc only
                self.fields = self.fields[:self.nnuc+1] 


                #percent cut for testing
                percent_test = 10
                N = len(react_data)

                # random seed
                randseed = 42

                Num_test  = int(N*percent_test/100)
                Num_train = N-Num_test

                train_set, test_set = torch.utils.data.random_split(react_data, [Num_train, Num_test], generator=torch.Generator().manual_seed(randseed))

                nbatch = 1 if device == torch.device('cpu') else torch.cuda.device_count()
                self.train_loader = DataLoader(dataset=train_set, batch_size=64*nbatch, shuffle=True)
                self.test_loader = DataLoader(dataset=test_set, batch_size=64*nbatch, shuffle=True)


    def train(self, model, optimizer, num_epochs, criterion, save_every_N=np.Inf):
            '''
            save_every_N - int representing every N epochs output the pytorch model.
                         Defaulted at infinity, meaning it won't output intermediately
            '''

            if save_every_N < np.Inf:
                os.mkdir(self.output_dir+'intermediate_output/')


            #As we test different loss functions, its important to keep a consistent one when
            #plotting or else we have no way of comparing them.
            criterion_plotting = nn.MSELoss()

            #plot storage
            self.cost_per_epoc = [] #stores total loss at each epoch
            self.component_losses_test = [] #stores component wise loss at each epoch (test data)
            self.component_losses_train = [] #stores component wise loss at each epoch (train data)
            self.cost_per_epoc_test = [] #stores total cost per epoc on testing data

            for epoch in range(num_epochs):
                losses = []
                plotting_losses = []

                for batch_idx, (data, targets) in enumerate(self.train_loader):
                    # Get data to cuda if possible
                    data = data.to(device=device)
                    targets = targets.to(device=device)

                    # forward
                    pred = model(data)
                    loss = criterion(pred, targets)

                    losses.append(loss.item())

                    if self.DO_PLOTTING:
                        with torch.no_grad():
                            loss_plot = criterion_plotting(pred, targets)
                            plotting_losses.append(loss_plot.item())

                            loss_c = component_loss_f(pred, targets)
                            loss_c = np.array(loss_c.tolist())
                            if batch_idx == 0:
                                component_loss = loss_c
                            else:
                                component_loss = component_loss + loss_c

                    # backward
                    optimizer.zero_grad()
                    loss.backward()

                    # gradient descent or adam step
                    optimizer.step()


                self.logger.write(f"Cost at epoch {epoch} is {sum(losses) / len(losses)}")
                self.component_losses_train.append(component_loss/batch_idx)
                self.cost_per_epoc.append(sum(plotting_losses) / len(plotting_losses))


                if self.DO_PLOTTING:
                    model.eval()

                    with torch.no_grad():
                        #Evaulate NN on testing data.
                        for batch_idx, (data, targets) in enumerate(self.test_loader):
                            # Get data to cuda if possible
                            data = data.to(device=device)
                            targets = targets.to(device=device)

                            # forward
                            pred = model(data)
                            loss = criterion_plotting(pred, targets)
                            losses.append(loss.item())

                            loss_c = component_loss_f(pred, targets)
                            loss_c = np.array(loss_c.tolist())
                            if batch_idx == 0:
                                component_loss = loss_c
                            else:
                                component_loss = component_loss + loss_c

                        self.cost_per_epoc_test.append(sum(losses) / len(losses))
                        self.component_losses_test.append(component_loss/batch_idx)

                    model.train()


                with torch.no_grad():
                    if epoch % save_every_N == 0 and epoch != 0:
                        directory = self.output_dir+'intermediate_output/epoch'+str(epoch)+'/'
                        os.mkdir(directory)

                        try:
                            state_dict = model.module.state_dict()
                        except AttributeError:
                            state_dict = model.state_dict()

                        torch.save(state_dict, directory+'my_model.pt')
                        np.savetxt(directory + "/cost_per_epoch.txt", self.cost_per_epoc)
                        np.savetxt(directory + "/component_losses_test.txt", self.component_losses_test)
                        np.savetxt(directory + "/component_losses_train.txt", self.component_losses_train)


                        plot_class = plotting_standard(model, self.fields, self.test_loader, self.cost_per_epoc, 
                                                       np.array(self.component_losses_test), np.array(self.component_losses_train), 
                                                       self.cost_per_epoc_test, directory)

                        plot_class.do_all_plots()


            self.component_losses_test = np.array(self.component_losses_test)
            self.component_losses_train = np.array(self.component_losses_train)



            self.model = model

            if self.SAVE_MODEL:
                self.logger.write("Saving...")
                file_name = self.output_dir + 'my_model.pt'
                if os.path.exists(file_name):
                    self.logger.write(f"Overwriting file: {file_name}")
                    os.rename(file_name, file_name+'.backup')

                try:
                    state_dict = model.module.state_dict()
                except AttributeError:
                    state_dict = model.state_dict()
                torch.save(state_dict, file_name)
                np.savetxt(self.output_dir + "/cost_per_epoch.txt", self.cost_per_epoc)
                np.savetxt(self.output_dir + "/component_losses_test.txt", self.component_losses_test)
                np.savetxt(self.output_dir + "/component_losses_train.txt", self.component_losses_train)


    def plot(self):
            self.logger.write("Plotting...")

            plot_class = plotting_standard(self.model, self.fields, self.test_loader, self.cost_per_epoc, self.component_losses_test,
                        self.component_losses_train, self.cost_per_epoc_test, self.output_dir, self.LOG_MODE)

            plot_class.do_all_plots()
