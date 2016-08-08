# New idea for comparisons: up until now, I have regridded all of the data
# onto a common grid in order to compute spatial correlations and biases.
# This always has the effect of smoothing the data. Perhaps a better way is
# to regrid by just selecting the points closest to the smallest grids points.
# I.e., instead of interpolating/averaging, at each target grid point, select
# the source grid point with coordinates closest to the target coordinates.
import numpy as numpy
import matplotlib as matplotlib
import matplotlib.pyplot as pyplot
import pylab as pylab
import mytaylor as taylor
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

inputdir = "~/cmor/climo"
varname = "cltisccp"

# we will need this to mask ice and land
maskprefix = "~/cmor/climo/misr/cltmisr/cltmisr_"

testinfo = numpy.array([
        ["cltmodis"     ,"modis"         ,"MODIS"          ,"DeepSkyBlue"],
        ["cltmodis_mask","modis"         ,"MODIS"          ,"RoyalBlue"  ],
        ["cltmisr"      ,"misr"          ,"MISR"           ,"DarkGreen"  ],
        ["cltisccp"     ,"isccp1984-1995","ISCCP 1984-1995","Firebrick"  ],
        ["cltisccp"     ,"isccp1996-2007","ISCCP 1996-2007","Crimson"    ],
    ])
cntlinfo = numpy.array([
        ["cltisccp","isccp","ISCCP"],
        ["cltisccp","isccp","ISCCP"],
        ["cltisccp","isccp","ISCCP"],
        ["cltisccp","isccp","ISCCP"],
        ["cltisccp","isccp","ISCCP"],
    ])
testvars = testinfo[:,0]
testnames = testinfo[:,1]
testlabels = testinfo[:,2]
testcolors = testinfo[:,3]
cntlvars = cntlinfo[:,0]
cntlnames = cntlinfo[:,1]
cntllabels = cntlinfo[:,2]
ntest = len(testnames)

# initialize statistics variables
ratio  = numpy.zeros([ntest])
cc     = numpy.zeros([ntest])
bias   = numpy.zeros([ntest])
rmserr = numpy.zeros([ntest])

# open figure
figure = pyplot.figure(figsize=(6.03,3.015))

# first, plot all points
for iplot in [0,1]:
    for i in range(ntest):

        # read data
        testdir = inputdir+"/"+testnames[i]+"/"+testvars[i]
        testprefix = testdir+"/"+testvars[i]+"_"
        test = myvars.taylor_var(testprefix,testvars[i])

        cntldir = inputdir+"/"+cntlnames[i]+"/"+cntlvars[i]
        cntlprefix = cntldir+"/"+cntlvars[i]+"_"
        cntl = myvars.taylor_var(cntlprefix,cntlvars[i])

        # regrid
        if len(test.lat) < len(cntl.lat) or len(test.lon) < len(cntl.lon):
            cntl.regrid(test.lat,test.lon)
        elif len(cntl.lat) < len(test.lat) or len(cntl.lon) < len(test.lon):
            test.regrid(cntl.lat,cntl.lon)

        # mask missing values
        test.mask_area(cntl)
        cntl.mask_area(test)

        # only ice-free ocean?
        if iplot == 1:
            maskvar = myvars.taylor_var(maskprefix,"cltmisr")
            maskvar.regrid(test.lat,test.lon)
            test.mask_area(maskvar)
            cntl.mask_area(maskvar)

        # subset region
        test.region_subset((-60.0,60.0),(0.0,360.0))
        cntl.region_subset((-60.0,60.0),(0.0,360.0))

        # calculate taylor statistics
        ratio[i],cc[i],bias[i],rmserr[i] = test.taylor_statistics(cntl)
        
    # make figure
    ax = figure.add_subplot(1,2,iplot+1,frameon=False)
    taylor_diagram = taylor.Taylor_diagram(
            ax,ratio,cc,bias,
            casecolors=testcolors,
        )

# test labels
figure.axes[0].text(
        1.17,0.30,"MODIS retrieval",
        color="DeepSkyBlue",
        horizontalalignment="right"
    )
figure.axes[0].text(
        1.02,0.60,"MODIS mask",
        color="RoyalBlue",
        horizontalalignment="right"
    )
figure.axes[0].text(
        0.83,0.47,"MISR",
        color="DarkGreen",
        horizontalalignment="right"
    )
figure.axes[0].text(
        0.95,0.12,r"ISCCP 1",
        color="Firebrick",
        horizontalalignment="right"
    )
figure.axes[0].text(
        1.05,0.11,r"ISCCP 2",
        color="Crimson",
        horizontalalignment="left"
    )

# plot labels
figure.axes[0].text(
        0.05,0.15,r"Total cloud",
        color="Black",
        horizontalalignment="left"
    )
figure.axes[1].text(
        0.05,0.15,r"Total cloud",
        color="Black",
        horizontalalignment="left"
    )
figure.axes[0].text(
        0.05,0.05,r"60S - 60N",
        color="Black",
        horizontalalignment="left",
    )
figure.axes[1].text(
        0.05,0.05,r"Ice-free ocean",
        color="Black",
        horizontalalignment="left",
    )

# adjust subplots
figure.subplots_adjust(left   = 0.01)
figure.subplots_adjust(bottom = 0.15)
figure.subplots_adjust(top    = 1.00)
figure.subplots_adjust(right  = 0.95)
figure.subplots_adjust(hspace = 0.05)
figure.subplots_adjust(wspace = 0.10)
 
# save figure
figure.savefig("../graphics/clt_retrievals_taylor.pdf",format="pdf")
