import yt
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import torch
import os

def make_movie(file_list, save_dir='./', var='enuc', movie_name="flame.mp4"):
    i = 1
    for file in file_list:
        ds = yt.load(file)
        sl = yt.SlicePlot(ds,2,var)
        sl.save("movie_imag{}.png".format(str(i).zfill(4)))
        i+=1
    os.system("ffmpeg -r 60 -pattern_type glob -i 'movie_imag*.png' -vcodec mpeg4 -y {}".format(movie_name))
    os.system("rm movie_imag*")
    Video("movie.mp4", embed=True)

class plotting_standard:
    #class to make it easy to plot things from driver. Set up all the data
    #then just call the methods for whatever plots you want.

    def __init__(self, model, fields, test_loader, cost_per_epoc,
                 component_losses_test, component_losses_train,
                 cost_per_epoc_test, output_dir, LOG_MODE=True):
        self.LOG_MODE = LOG_MODE

        self.model = model
        self.fields = fields
        self.nnuc = len(fields)-1
        self.N_fields = len(fields)
        self.test_loader = test_loader
        self.cost_per_epoc = cost_per_epoc
        self.component_losses_test = component_losses_test
        self.component_losses_train = component_losses_train
        self.cost_per_epoc_test = cost_per_epoc_test
        self.output_dir = output_dir

        isdir = os.path.isdir(output_dir)
        if not isdir:
            os.mkdir(output_dir)


    def do_prediction_vs_solution_plot(self):
        ############ Prediction Vs Solution Plot Should fall one y=x line.

        plt.figure()
        #N = react_data.output_data.shape[1]
        colors = matplotlib.cm.rainbow(np.linspace(0, 1, len(self.fields)))
        #fields = [field[1] for field in yt.load(react_data.output_files[0]).field_list]
        self.model.eval()

        with torch.no_grad():
            losses = []

            for batch_idx, (data, targets) in enumerate(self.test_loader):
            #pulling this data out and storing it then calling matplotlib as
            #few times as possible is much faster.
                if batch_idx==0:
                    data_whole = data
                    targets_whole = targets
                else:
                    data_whole = torch.cat((data_whole, data))
                    targets_whole = torch.cat((targets_whole, targets))

            data_whole = data_whole.cuda()
                    
            #for batch_idx, (data, targets) in enumerate(self.test_loader):
            pred = self.model(data_whole)

            if self.LOG_MODE:
                # convert all mass fractions back from their log form
                data_whole[:,:self.nnuc] = torch.exp(-1.0/data_whole[:,:self.nnuc])
                targets_whole[:,:self.nnuc] = torch.exp(-1.0/targets_whole[:,:self.nnuc])
                pred[:,:self.nnuc] = torch.exp(-1.0/pred[:,:self.nnuc])
            
            pred = pred.cpu()
            targets_whole = targets_whole.cpu()
                
            colors_nm1 = np.tile(colors, (pred.shape[0]-1, 1))
            for j in range(pred.shape[1]):
                plt.scatter(pred[0,j], targets_whole[0,j], color=colors[j], label=self.fields[j])
            plt.scatter(pred[1:, :], targets_whole[1:, :], c=colors_nm1)
            
            plt.plot(np.linspace(0, 1), np.linspace(0,1), '--', color='orange')
            #plt.legend(yt.load(react_data.output_files[0]).field_list, colors=colors)
            plt.legend(bbox_to_anchor=(1, 1))
            plt.xlabel('Prediction')
            plt.ylabel('Solution')
            plt.savefig(self.output_dir + "/prediction_vs_solution.png", bbox_inches='tight')

            plt.yscale("log")
            plt.xscale("log")
            plt.ylim([1.e-16, 1.e1])
            plt.xlim([1.e-16, 1.e3])
            plt.savefig(self.output_dir + "/prediction_vs_solution_log.png", bbox_inches='tight')

        self.model.train()

        plt.close()


    def do_cost_per_epoch_plot(self):
        ############## Cost per eppoc plot#####################
        fig = plt.figure()
        gs = fig.add_gridspec(2,1, hspace=0)
        axs = gs.subplots(sharex=True)
        epocs = np.linspace(1, len(self.cost_per_epoc), num=len(self.cost_per_epoc))

        axs[0].plot(epocs, self.cost_per_epoc, label='Training Data')
        axs[0].plot(epocs, self.cost_per_epoc_test, label='Testing Data')

        axs[1].semilogy(epocs, self.cost_per_epoc)
        axs[1].semilogy(epocs, self.cost_per_epoc_test)

        fig.suptitle('Overall cost of training data')
        axs[1].set_xlabel("Num Epochs")
        axs[1].set_ylabel('Log Cost')
        axs[0].set_ylabel('Cost (MSE)')
        # Hide x labels and tick labels for all but bottom plot.
        for ax in axs:
            ax.label_outer()
        axs[0].legend(bbox_to_anchor=(1, 1))
        fig.savefig(self.output_dir + "/cost_vs_epoch.png", bbox_inches='tight')

        plt.close(fig)


    def do_component_loss_train_plot(self):

        #Component losses  train
        fig = plt.figure()
        gs = fig.add_gridspec(2,1, hspace=0)
        axs = gs.subplots(sharex=True)

        N = self.component_losses_train.shape[0]
        for i in range(self.component_losses_train.shape[1]):
            axs[0].plot(np.linspace(1, N, num=N),
                        self.component_losses_train[:, i], label=self.fields[i])

        for i in range(self.component_losses_train.shape[1]):
            axs[1].semilogy(np.linspace(1, N, num=N),
                            self.component_losses_train[:, i], label=self.fields[i])

        fig.suptitle('Component wise error in training data')
        axs[1].set_xlabel("Num Epochs")
        axs[1].set_ylabel('Log Cost')
        axs[0].set_ylabel('Cost (MSE)')
        # Hide x labels and tick labels for all but bottom plot.
        for ax in axs:
            ax.label_outer()
        plt.legend(bbox_to_anchor=(1, 2))
        fig.savefig(self.output_dir + "/component_training_loss.png", bbox_inches='tight')

        plt.close(fig)


    def do_component_loss_test_plot(self):
        #Component losses  test
        fig = plt.figure()
        gs = fig.add_gridspec(2,1, hspace=0)
        axs = gs.subplots(sharex=True)

        N = self.component_losses_test.shape[0]
        for i in range(self.component_losses_test.shape[1]):
            axs[0].plot(np.linspace(1, N, num=N), self.component_losses_test[:, i],
                        label=self.fields[i])

        for i in range(self.component_losses_test.shape[1]):
            axs[1].semilogy(np.linspace(1, N, num=N), self.component_losses_test[:, i],
                            label=self.fields[i])

        fig.suptitle('Component wise error in testing data')
        axs[1].set_xlabel("Num Epochs")
        axs[1].set_ylabel('Log Cost')
        axs[0].set_ylabel('Cost (MSE)')
        # Hide x labels and tick labels for all but bottom plot.
        for ax in axs:
            ax.label_outer()
        plt.legend(bbox_to_anchor=(1, 2))
        fig.savefig(self.output_dir + "/component_testing_loss.png", bbox_inches='tight')

        plt.close(fig)


    def do_all_plots(self):
        self.do_cost_per_epoch_plot()
        self.do_component_loss_train_plot()
        self.do_component_loss_test_plot()
        self.do_prediction_vs_solution_plot()
