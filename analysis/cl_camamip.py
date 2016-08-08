import numpy as numpy
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
 
varname = "cl"
testinfo = numpy.array([
        ["cam3_amip","CAM3 T42","SeaGreen"],
        ["cam3_amip_fv","CAM3 FV","MediumAquamarine"],
        ["cam4_1deg_release_amip","CAM4","RoyalBlue"]
    ])
testnames = testinfo[:,0]
testlabels = testinfo[:,1]
testcolors = testinfo[:,2]
ntest = len(testnames)

# open figure
golden = (1+numpy.sqrt(5.0))/2.0
width = 6.03
height = width/golden
figure = pyplot.figure(figsize=(width,height))
ax = figure.add_subplot(111)

for i in range(ntest):

    # read data
    directory = inputdir+"/"+testnames[i]+"/"+varname
    filename = directory+"/"+varname+"_ANN.nc"
    test = myvars.cl_var(filename,varname)

    # calculate global average
    cl = test.weighted_area_average()

    # plot control data
    ax.plot(
            cl,test.lev,'.',
            color=testcolors[i],
            linestyle="solid",
            linewidth=2.0,
            markersize=10.0,
            label=testlabels[i]
        )

# set axes ticks and labels
ax.set_ylim(ax.get_ylim()[::-1])
ax.set_xlabel("Cloud amount (\%)")
ax.set_ylabel("Hybrid level")
ax.grid()

# add a legend
pyplot.legend(loc="lower left",handlelength=3)

# adjust subplots
figure.subplots_adjust(left=0.1)
figure.subplots_adjust(bottom=0.15)
figure.subplots_adjust(top=0.95)
figure.subplots_adjust(right=0.95)

# save figure
figure.savefig("../graphics/cl_camamip.pdf",format="pdf")
