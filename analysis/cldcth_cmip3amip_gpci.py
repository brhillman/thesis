import numpy as numpy
import Scientific.IO.NetCDF as netcdf
import matplotlib.pyplot as pyplot
import pylab as pylab
import myvars as myvars

# fix fonts
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

# case and variable parameters
inputdir = "~/cmor/climo"
maskfile = inputdir+"/misr/clMISR/clMISR_JJA.nc"
maskvar = myvars.clcth_var(maskfile,"clMISR")
 
testinfo = numpy.array([
        ["misr","MISR"],
        ["cam3_amip","CAM3"],
        ["am2_amip","AM2" ],
    ])
varname = "clMISR"
cth_bnds = numpy.array([
        0.0,0.5,1.0,1.5,2.0,2.5,3.0,4.0,5.0,
        7.0,9.0,11.0,13.0,15.0,17.0,23.0
    ])
testnames = testinfo[:,0]
testlabels = testinfo[:,1]
ntest = len(testnames)

# figure size
golden = (1+numpy.sqrt(5.0))/2.0
width = 6.03
height = width/golden
figure = pyplot.figure(figsize=(width,height))
cmap = pyplot.cm.get_cmap("Blues",20)

# cross section specification
lats = numpy.arange(-1.0,35.1,3)
lons = numpy.arange(187.0,235.1,4)
#lats = numpy.arange(35.0,-1.1,-3)
#lons = numpy.arange(235.0,186.9,-4)

# loop over test cases
for itest in range(ntest):

    # read data
    directory = inputdir+"/"+testnames[itest]+"/"+varname
    filename = directory+"/"+varname+"_ANN.nc"
    test = myvars.clcth_var(filename,varname)
    test.calculate()

    # grab points along cross-section
    lat = test.lat
    lon = test.lon
    histogram = numpy.zeros([len(test.cth),len(lats)])
    for idx in range(len(lats)):
        latind = numpy.abs(lat-lats[idx]).argmin()
        lonind = numpy.abs(lon-lons[idx]).argmin()
        histogram[:,idx] = test.histogram[:,latind,lonind]

    # make plot
    ax = figure.add_subplot(2,2,itest+1)
    
    plot = ax.pcolor(
            histogram,
            vmin=0.0,vmax=30.0,
            cmap=cmap,
        )
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    ax.text(
            12.5,14.0,testlabels[itest],
            horizontalalignment="right",
            verticalalignment="top"
        )

# ticks and tick labels
cthticklabels = []
cthticks = numpy.arange(len(cth_bnds))
for i in range(len(cth_bnds)):
    cthticklabels.append("%.1f"%cth_bnds[i])

latticklabels = []
latticks = numpy.arange(len(lats))+0.5
for i in range(len(lats)):
    latticklabels.append("%.1f"%lats[i])

# now adjust ticks
for i in range(ntest):
    ax = figure.add_axes(figure.axes[i])
    ax.set_yticks(cthticks[::2])
    ax.set_xticks(latticks[::3])
    ax.set_ylim([0,len(cth_bnds)-1])
    ax.set_xlim([0,len(lats)])

for i in [0,2]:
    ax = figure.add_axes(figure.axes[i])
    ax.set_yticklabels(cthticklabels[::2])
    ax.set_ylabel("CTH (km)")

for i in [1,2]:
    ax = figure.add_axes(figure.axes[i])
    ax.set_xticklabels(latticklabels[::3])
    ax.set_xlabel("Latitude")

cax = figure.add_axes([0.555,0.30,0.435,0.04])
cb = pyplot.colorbar(
        plot,
        orientation='horizontal',
        cax=cax,
        drawedges=False,
    )
cb.set_label(r"Frequency (\%)")


# adjust
figure.subplots_adjust(left=0.1)
figure.subplots_adjust(bottom=0.10)
figure.subplots_adjust(top=0.99)
figure.subplots_adjust(right=0.99)
figure.subplots_adjust(wspace=0.05)
figure.subplots_adjust(hspace=0.05)

# save figure
figure.savefig("../graphics/cldcth_cmip3amip_gpci.pdf",format="pdf")
