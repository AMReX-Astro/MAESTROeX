import sys
import yt
yt.funcs.mylog.setLevel(50)

def plot_file(filename, var):

    ds = yt.load(filename)
    fig = yt.SlicePlot(ds, 'z', var, center=[ds.domain_center[0],ds.domain_center[1],0])
    fig.set_log(var, True)
    plotname = filename.split('/')[-1] + "_{}.png".format(var)
    if var == "Hnuc":
        fig.set_zlim(var, zmin=1.e6, zmax=1.e15)
    print(plotname)

    fig.save(plotname)

if __name__ == "__main__":
    if len(sys.argv) > 2:
        filename, var = sys.argv[1:3]
    else:
        filename = sys.argv[1]
        var = "magvel"

    plot_file(filename, var)
