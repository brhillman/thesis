import numpy as numpy
import matplotlib as matplotlib
import matplotlib.pyplot as pyplot
from mpl_toolkits.basemap import Basemap
import pylab as pylab
import myvars as myvars

params = {
        'axes.titlesize': 10,
        'axes.labelsize': 10,
        'text.fontsize': 10,
        'legend.fontsize': 10,
        'xtick.labelsize': 10,
        'ytick.labelsize': 10,
        'font.family': 'serif',
        'font.family': 'serif',
        'font.size': 10,
        'text.usetex': True
    }
pylab.rcParams.update(params)

input_dir = "~/cmor/climo"
varname = "swcre"

# test and control case parameters
testinfo = numpy.array([
        ["cam3_hot4k","CAM3 HOT4K"],
        ["am2_hot4k" ,"AM2 HOT4K" ],
    ])
cntlinfo = numpy.array([
        ["cam3_amip","CAM3 AMIP"],
        ["am2_amip" ,"AM2 AMIP" ]
    ])
testnames = testinfo[:,0]
testlabels = testinfo[:,1]
cntlnames = cntlinfo[:,0]
cntllabels = cntlinfo[:,1]
ntest = len(testnames)

# open figure
figure = pyplot.figure(figsize=(6.03,2.35))

# contour levels
cmap = pyplot.cm.RdBu_r
clevels = numpy.arange(-30,30.1,2.5)

for i in range(ntest):

    # read data
    testdir = input_dir+"/"+testnames[i]+"/"+varname
    testfile = testdir+"/"+varname+"_ANN.nc"
    test = myvars.map_var(filename=testfile,varname=varname)

    cntldir = input_dir+"/"+cntlnames[i]+"/"+varname
    cntlfile = cntldir+"/"+varname+"_ANN.nc"
    cntl = myvars.map_var(filename=cntlfile,varname=varname)

    # regrid
    if len(test.lat) < len(cntl.lat) or len(test.lon) < len(cntl.lon):
        cntl.regrid(test.lat,test.lon)
    elif len(cntl.lat) < len(test.lat) or len(cntl.lon) < len(test.lon):
        test.regrid(cntl.lat,cntl.lon)

    # mask missing values
    test.mask_area(cntl)
    cntl.mask_area(test)

    # calculate bias
    bias = test.copy()
    bias.data = test.data - cntl.data

    # calculate means and standard deviation
    bias_mean = bias.weighted_area_average()
    bias_std = bias.weighted_std()

    # make subplot
    ax = figure.add_subplot(1,2,i+1)

    # transform coordinates
    datamap = Basemap(projection='robin',lon_0=0.5*(bias.lon[0]+bias.lon[-1]))
    x,y = datamap(*numpy.meshgrid(bias.lon,bias.lat))

    # plot data
    plot  = datamap.contourf(x,y,bias.data,clevels,extend="both",cmap=cmap)
    pyplot.title(testlabels[i]+" - "+cntllabels[i]+" (%.0f %s)"%(bias_mean,r"${\rm W}/{\rm m}^2$"))

    # draw lines
    datamap.drawcoastlines(linewidth=0.5)
    datamap.drawmapboundary(linewidth=0.5)

# adjust subplots
figure.subplots_adjust(left   = 0.00)
figure.subplots_adjust(bottom = 0.15)
figure.subplots_adjust(top    = 1.00)
figure.subplots_adjust(right  = 1.00)
figure.subplots_adjust(hspace = 0.05)
figure.subplots_adjust(wspace = 0.01)

# common colorbar
ticks = [-20,-10,0,10,20]
cax = figure.add_axes([0.10,0.15,0.8,0.075])
cb = pyplot.colorbar(plot,
                     orientation='horizontal',cax=cax,
                     drawedges=False,ticks=ticks)
cb.set_label(r"Cloud radiative effect (${\rm W}/{\rm m}^2$)")

# save figure
figure.savefig("../graphics/swcrechanges_cmip3hot4k_map.pdf",format="pdf")
