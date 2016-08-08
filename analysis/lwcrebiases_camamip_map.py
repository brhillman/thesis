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

basedir = "~/cmor/climo"
varname = "lwcre"

# test and control case parameters
test_info = numpy.array([
        ["cam3_amip"             ,"CAM3 T42"],
        ["cam3_amip_fv"          ,"CAM3 FV" ],
        ["cam4_1deg_release_amip","CAM4"],
        ["cam5_1deg_release_amip","CAM5"],
    ])
cntl_info = numpy.array([
        ["ceres-amwg"            ,"CERES"],
        ["ceres-amwg"            ,"CERES"],
        ["ceres-amwg"            ,"CERES"],
        ["ceres-amwg"            ,"CERES"],
    ])
test_names = test_info[:,0]
test_labels = test_info[:,1]
cntl_names = cntl_info[:,0]
cntl_labels = cntl_info[:,1]
ntest = len(test_names)

# open figure
figure = pyplot.figure(figsize=(6.03,4.25))

# contour levels
cmap = pyplot.cm.RdBu_r
clevels = numpy.arange(-30,30.1,2.5)

for i in range(ntest):

    # read data
    testdir = basedir+"/"+test_names[i]+"/"+varname
    testfile = testdir+"/"+varname+"_ANN.nc"
    test = myvars.map_var(filename=testfile,varname=varname)

    cntldir = basedir+"/"+cntl_names[i]+"/"+varname
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
    ax = figure.add_subplot(2,2,i+1)

    # transform coordinates
    datamap = Basemap(projection='robin',lon_0=0.5*(bias.lon[0]+bias.lon[-1]))
    x,y = datamap(*numpy.meshgrid(bias.lon,bias.lat))

    # plot data
    plot  = datamap.contourf(x,y,bias.data,clevels,extend="both",cmap=cmap)
    pyplot.title(test_labels[i]+" - "+cntl_labels[i]+" (%.0f %s)"%(bias_mean,r"${\rm W}/{\rm m}^2$"))

    # draw lines
    datamap.drawcoastlines(linewidth=0.5)
    datamap.drawmapboundary(linewidth=0.5)

# adjust subplots
figure.subplots_adjust(left   = 0.00)
figure.subplots_adjust(bottom = 0.12)
figure.subplots_adjust(top    = 0.98)
figure.subplots_adjust(right  = 1.00)
figure.subplots_adjust(hspace = 0.02)
figure.subplots_adjust(wspace = 0.01)

# common colorbar
ticks = [-20,-10,0,10,20]
cax = figure.add_axes([0.10,0.085,0.8,0.045])
cb = pyplot.colorbar(plot,
                     orientation='horizontal',cax=cax,
                     drawedges=False,ticks=ticks)
cb.set_label(r"Cloud radiative effect (${\rm W}/{\rm m}^2$)")

# save figure
figure.savefig("../graphics/lwcrebiases_camamip_map.pdf",format="pdf")
