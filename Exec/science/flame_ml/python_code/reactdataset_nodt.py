
from torch.utils.data.dataset import Dataset
import yt
import numpy as np
import matplotlib.pyplot as plt
import os
from IPython.display import Video
from glob import glob
import torch
import warnings
import sys

import pandas as pd
from torch.utils.data.dataset import Dataset
from torch.utils.data import DataLoader
yt.funcs.mylog.setLevel(40) # Gets rid of all of the yt info text, only errors.
import matplotlib.cbook
warnings.filterwarnings("ignore",category=matplotlib.cbook.mplDeprecation) #ignore plt depreciations
import gc

class ReactDataset(Dataset):

    def __init__(self, data_path, input_prefix, output_prefix, plotfile_prefix, DEBUG_MODE=False, data_range=3):
        #loading data
        #Load input and output data
        self.DEBUG_MODE = DEBUG_MODE
        self.input_prefix = input_prefix
        self.output_prefix = output_prefix
        self.data_path = data_path

        self.do_flame_cut = True
        self.xbeg = data_range * 0.1
        self.xend = self.xbeg + 0.01
        self.nf = 2   # how many extra local flame data to append to dataset (<= 6)
        if self.nf > 0:
            self.xrange = [i * 0.1 for i in range(1,7)]

        self.input_files  = self.get_files(data_path, input_prefix)
        self.output_files = self.get_files(data_path, output_prefix)

        print("Loading Input Files...")
        self.input_data, self.input_break_file  = self.load_files(self.input_files, inputs=True)
        print("Loading Output Files...")
        self.output_data, self.output_break_file = self.load_files(self.output_files)

        #we want them to be the same length - so cut off data if we need to

        if (self.input_break_file is None) or (self.output_break_file is None):
            ind_in = len(self.input_files)
            ind_out = len(self.output_files)
        else:
            ind_in = self.input_files.index(self.input_break_file)
            ind_out = self.output_files.index(self.output_break_file)


        if ind_in <= ind_out:
            ind = ind_in
        else:
            ind = ind_out

        #cut excess.
        self.input_data = self.input_data[0:ind, :, :]
        self.output_data = self.output_data[0:ind, :, :]
        self.input_files = self.input_files[0:ind]
        self.output_files = self.output_files[0:ind]

        #get max enuc
        self.n_output = int(self.output_data.shape[1])//2
        self.enuc_fac = torch.max(self.output_data[:, self.n_output-1, :])

        print("Loaded data successfully!")

    def get_files(self, data_path, prefix):
        data_files = glob(data_path + prefix)
        data_files = sorted(data_files)
        for data in data_files:
            if data[-7:] == 'endstep':
                data_files.remove(data)

        if self.DEBUG_MODE:
            data_files = data_files[:5]
        return data_files


    def load_files(self, file_list, inputs=False):

        break_file = None

        #Store data each row corresponds to data acros the grid of a different field.
        for j, file in enumerate(file_list):

            try:
                ds = yt.load(file)
                dt = ds.current_time.to_value()
                #Store data each row corresponds to data acros the grid of a different field.
                if self.do_flame_cut:
                    flame_loc = self.get_flame_loc(file)
                    half = yt.YTArray(0.5, 'cm')
                    lbound = yt.YTArray(0.1, 'cm')
                    if flame_loc-half > lbound:
                        lbound = flame_loc-half
                    ad_flame = ds.r[self.xbeg:self.xend, lbound:flame_loc+half]
                    if self.nf > 0:
                        ad_flame = ad_flame + ds.r[self.xrange[0]:self.xrange[0] + 0.01, lbound:flame_loc+half]
                        for n in range(1,self.nf):
                            ad_flame = ad_flame + ds.r[self.xrange[n]:self.xrange[n] + 0.01, lbound:flame_loc+half]
                    ad = ds.r[self.xbeg:self.xend, :]
                else:
                    ad_flame = []
                    ad = ds.all_data()

                for i,field in enumerate(ds._field_list):
                    # add repeating data of the flame itself (nf times)
                    if i == 0:
                        data = np.zeros([len(ds._field_list), len(ad[field]) + len(ad_flame[field])])

                    if self.nf > 0:
                        data[i,:] = np.concatenate((np.array(ad[field]), np.array(ad_flame[field])))
                    else:
                        data[i,:] = np.array(ad[field])
            except:
                pass


            if j == 0:
                #dt
                if inputs:
                    dt_tensor = dt*torch.ones([1,1,data.shape[1]])
                data = torch.from_numpy(data.reshape((1,data.shape[0],data.shape[1])))
                if inputs:
                    data = torch.cat((dt_tensor,data), dim=1)

                data_set = data
            else:

                try:
                    #dt
                    NUM_GRID_CELLS = data_set.shape[2]
                    if inputs:
                        dt_tensor = dt*torch.ones([1,1,data.shape[1]])

                    data = torch.from_numpy(data.reshape((1,data.shape[0],data.shape[1])))
                    #print(data.shape)
                    if inputs:
                        data = torch.cat((dt_tensor,data), dim=1)

                    #If we have more data - cut data
                    if data.shape[2] > NUM_GRID_CELLS:
                        data = data[:, : , :NUM_GRID_CELLS]
                    #We need to get more data.
                    elif data.shape[2] < NUM_GRID_CELLS:
                        #double size of cut
                        if flame_loc-2*half > lbound:
                            lbound = flame_loc-2*half
                        ad_flame = ds.r[self.xbeg:self.xend, lbound:flame_loc+2*half]
                        if self.nf > 0:
                            ad_flame = ad_flame + ds.r[self.xrange[0]:self.xrange[0] + 0.01, lbound:flame_loc+2*half]
                            for n in range(1,self.nf):
                                ad_flame = ad_flame + ds.r[self.xrange[n]:self.xrange[n]+0.01, lbound:flame_loc+2*half]
                        ad = ds.r[self.xbeg:self.xend, :]
                        for i,field in enumerate(ds._field_list):
                            if i == 0:
                                data = np.zeros([len(ds._field_list), len(ad[field]) + len(ad_flame[field])])

                            if self.nf > 0:
                                data[i,:] = np.concatenate((np.array(ad[field]), np.array(ad_flame[field])))
                            else:
                                data[i,:] = np.array(ad[field])
                        data = torch.from_numpy(data.reshape((1,data.shape[0],data.shape[1])))
                        if inputs:
                            dt_tensor = dt*torch.ones([1,1,NUM_GRID_CELLS])
                        data = data[:, : , :NUM_GRID_CELLS]
                        if inputs:
                            data = torch.cat((dt_tensor,data), dim=1)


                    #z2 = data.reshape((1,data.shape[0],data.shape[1]))
                    #print(data.size())
                    data_set = torch.cat((data_set, data))
                    #torch.torch.stack([data_set,data], dim=0)
                except:
                    print('invalid data in file stopping here: ', file)
                    break_file = file
                    break

        return data_set, break_file


    def __getitem__(self, index):
        #indexing data dataset[0]

        #indexing goes across domain first, then to next data set.

        file_number = int(np.floor(index/self.input_data.shape[2]))
        cell_number = int(index%self.input_data.shape[2])
        iout = [i for i in range(self.n_output)]
#         if self.n_output == 4:
#             iout = [0, 2, 3]

        #data
        X = self.input_data[file_number, 1:, cell_number] #no dt
        #labels
        Y = self.output_data[file_number, iout, cell_number]

        return (X.float(),Y.float())

    def __len__(self):
        #Each cell contains a new training value.
        #The data is of the form [N1,N2,N3]
        #where N1 is the number of plot files we loaded
        #N2 is the inputs we will use in the model (state / thermo info)
        #N3 is the number of cells in that plotfile. Note this is 2D and we might
        #not want to store this as a 1d array in the future.

        #But for now we're simply just trying to learn a mapping between the input
        #thermo state and the output thermo state.
        return self.input_data.shape[2]* self.input_data.shape[0]


    def cut_data_set(self,N):
        #We have about 8gb of data and i can't train with that much when we're just testing.
        self.input_data = self.input_data[1:N,:,:]
        self.output_data = self.output_data[1:N,:,:]


    def get_flame_loc(self, file):
        if self.output_prefix[:-1] in file:
            output_file = True
            ds = yt.load(file)

        elif self.input_prefix[:-1] in file:
            #if we're in an input file, we must load corresponding output file to get
            #max hnuc and take that same clipping from the domain.
            output_file = False
            corr_out_file = self.data_path+self.output_prefix[:-1]+file[-9:]
            ds = yt.load(corr_out_file)
        else:
            print("Error, file must match eiter input or output prefix")
            sys.exit()

        #Cutting data
        #ds = yt.load(plotfiles[4])
        left  = ds.domain_left_edge
        right = ds.domain_right_edge



        ad = ds.all_data()
        x,y,z = ad.argmax("enuc")
        #Hnucs.append(ad.quantities.extrema("Hnuc"))
        # xs.append(x)
        # ys.append(y)
        # zs.append(z)


        #Another option is to average the ys along the x and take the max of that.
        # To do this we create bins along the x direction
        #The following bins by y, and with a weighted field of None,
        #For each y, we are just summing up all of enuc horizontally across.
        profile = ad.profile(("index", "y"), "enuc", weight_field=None, n_bins = 1000)
        #index of the max value of Hnuc
        index = profile.items()[0][1].argmax()
        #location of the max value
        flame_loc = profile.x[index]
        # rs.append(flame_loc)

        half = yt.YTArray(0.5, 'cm')
        lbound = yt.YTArray(0.1, 'cm')
        if flame_loc-half > lbound:
            lbound = flame_loc-half
        flame_slice = ds.r[self.xbeg:self.xend, lbound:flame_loc+half]
        #ds.domain_width.value
        #plot with slice this still plots the whole domain but zeros out everything else...

        #yt.SlicePlot(ds, 'z', 'rho', data_source = flame_slice, axes_unit='cm')


        #then we must take this slice and apply it to our inputfile
        if not output_file:
            ds = yt.load(file)
            flame_slice = ds.r[self.xbeg:self.xend, lbound:flame_loc+half]

        return flame_loc
