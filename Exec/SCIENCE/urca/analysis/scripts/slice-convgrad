#!/usr/bin/env python
import yt
from yt import derived_field
import numpy as np
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('infile', type=str, help='Name of input plotfile.')
parser.add_argument('-adiabatic', '--adiabatic', action='store_true', help='If supplied, plot adiabatic excess.')
parser.add_argument('-ledoux', '--ledoux', action='store_true', help='If supplied, plot ledoux excess.')
parser.add_argument('-convtype', '--convtype', action='store_true', help='If supplied, plot convective type.')
parser.add_argument('-w', '--width', type=float,
                    help='Width of slice (cm). Default is domain width.')
parser.add_argument('-sign', '--sign', action='store_true', help='If supplied, plot only the sign of the excess.')
parser.add_argument('-fmin', '--field_min', type=float, help='Minimum scale for colorbar.')
parser.add_argument('-fmax', '--field_max', type=float, help='Maximum scale for colorbar.')
parser.add_argument('-log', '--logscale', action='store_true', help='If supplied, use a log scale for the field.')
parser.add_argument('-symlog', '--symlog', action='store_true', help='If supplied, use symlog scaling, which is linear near zero, to accomodate positive and negative values.')
parser.add_argument('-linthresh', '--linthresh', type=float, help='Linear threshold for symlog scaling. (Default is 0.1)')
parser.add_argument('-dc', '--drawcells', action='store_true', help='If supplied, draw the cell edges.')
parser.add_argument('-dg', '--drawgrids', action='store_true', help='If supplied, draw the grids.')
parser.add_argument('-octant', '--octant', action='store_true', help='Sets slice view appropriately for octant dataset.')
args = parser.parse_args()

@derived_field(name='adiabatic_excess')
def _adiabatic_excess(field, data):
    return data[('boxlib','actual gradient')] - data[('boxlib','adiabatic gradient')]

@derived_field(name='ledoux_excess')
def _ledoux_excess(field, data):
    return data[('boxlib','actual gradient')] - data[('boxlib','ledoux gradient')]

@derived_field(name='sign_adiabatic_excess')
def _sign_adiabatic_excess(field, data):
    return np.sign(data[('boxlib','actual gradient')] - data[('boxlib','adiabatic gradient')])

@derived_field(name='sign_ledoux_excess')
def _sign_ledoux_excess(field, data):
    return np.sign(data[('boxlib','actual gradient')] - data[('boxlib','ledoux gradient')])

# Convective type is defined as follows:
# 1.0 = Convectively unstable wrt Ledoux and Adiabatic
# 0.0 = Semiconvection, stable wrt Ledoux but unstable wrt Adiabatic or v.v.
# -1.0 = Convectively stable wrt Ledoux and Adiabatic
@derived_field(name='conv_type')
def _conv_type(field, data):
    return 0.5*data['sign_ledoux_excess']+0.5*data['sign_adiabatic_excess']


def doit(field):
    ds = yt.load(args.infile)

    if not args.width:
        width = max(ds.domain_width)
    else:
        width = yt.YTQuantity(args.width, 'cm')

    maxv = ds.all_data().max(field)
    minv = ds.all_data().min(field)
    pos_maxv = np.ceil(np.log10(maxv))
    neg_maxv = np.ceil(np.log10(minv))
    logmaxv = max(pos_maxv, neg_maxv)
    linminv = min(abs(maxv), abs(minv))

    if args.octant:
        dcenter = width.in_units('cm').v/2.0
        cpos    = ds.arr([dcenter, dcenter, dcenter], 'cm')
        s = yt.SlicePlot(ds, 'x', field, center=cpos, width=width, origin="native")
    else:
        s = yt.SlicePlot(ds, 'x', field, center='c', width=width)

    s.annotate_scale()

    if args.drawcells:
        s.annotate_cell_edges()

    if args.drawgrids:
        s.annotate_grids()

    s.set_buff_size(2048)
    
    if minv < 0.0 and maxv > 0.0:
        s.set_cmap(field, 'PiYG')
        linthresh = 0.1
        if ((args.logscale or args.symlog) and
            not args.sign and not field=='conv_type'):
            if args.linthresh:
                linthresh = args.linthresh
            s.set_log(field, args.logscale, linthresh=linthresh)
        else:
            s.set_log(field, args.logscale)
        if args.sign or field=='conv_type':
            s.set_zlim(field, -1.0, 1.0)
            plot = s.plots[field]
            colorbar = plot.cb
            s._setup_plots()
            if field != 'conv_type':
                colorbar.set_ticks([-1, 0, 1])                
                colorbar.set_ticklabels(['$-1$', '$0$', '$+1$'])
            else:
                colorbar.set_ticks([-0.8, 0, 0.65])
                colorbar.ax.tick_params(axis=u'both', which=u'both',length=0)
                colorbar.ax.set_yticklabels(['stable', 'semiconvective', 'convective'], rotation=90)
        elif args.field_min and args.field_max:
            s.set_zlim(field, args.field_min, args.field_max)
        elif args.field_min and args.field_min < 0.0:
            s.set_zlim(field, args.field_min, -args.field_min)
        elif args.field_max and args.field_max > 0.0:
            s.set_zlim(field, -args.field_max, args.field_max)
        else:
            s.set_zlim(field, -linthresh, linthresh)
    else:
        s.set_cmap(field, 'Greens')    

    s.save('{}.slice.{}.png'.format(args.infile, field))

if __name__=='__main__':
    if args.convtype:
        field = 'conv_type'
        doit(field)
    if args.adiabatic:
        field = 'adiabatic_excess'
        if args.sign:
            field = 'sign_' + field
        doit(field)
    if args.ledoux:
        field = 'ledoux_excess'
        if args.sign:
            field = 'sign_' + field
        doit(field)
    if not args.adiabatic and not args.ledoux and not args.convtype:
        if args.sign:
            doit('sign_adiabatic_excess')
            doit('sign_ledoux_excess')
        else:
            doit('adiabatic_excess')
            doit('ledoux_excess')
        doit('conv_type')
