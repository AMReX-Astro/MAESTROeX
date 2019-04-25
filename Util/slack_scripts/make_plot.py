import sys
import yt
yt.funcs.mylog.setLevel(50)

def plot_file(filename, var):

    ds = yt.load(filename)
    fig = yt.SlicePlot(ds, 'z', var)
    plotname = filename.split('/')[-1] + ".png"

    print(plotname)

    fig.save(plotname)

if __name__ == "__main__":
    if len(sys.argv) > 2:
        filename, var = sys.argv[1:3]
    else:
        filename = sys.argv[1]
        var = "magvel"

    plot_file(filename, var)
