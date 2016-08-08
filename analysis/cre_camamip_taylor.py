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
varinfo = numpy.array([
        ["swcre","ceres-amwg","CERES-EBAF"],
        ["lwcre","ceres-amwg","CERES-EBAF"],
    ])
testinfo = numpy.array([
        ["erbe-amwg","ERBE","0.50"],
        ["cam3_amip","CAM3","SeaGreen"],
        ["cam4_1deg_release_amip","CAM4","RoyalBlue"],
        ["cam5_1deg_release_amip","CAM5","Firebrick"],
    ])
varnames = varinfo[:,0]
cntlnames = varinfo[:,1]
cntllabels = varinfo[:,2]
nvars = len(varnames)

testnames = testinfo[:,0]
testlabels = testinfo[:,1]
testcolors = testinfo[:,2]
ntest = len(testnames)

# initialize statistics variables
ratio  = numpy.zeros([ntest])
cc     = numpy.zeros([ntest])
bias   = numpy.zeros([ntest])
rmserr = numpy.zeros([ntest])

# open figure
figure = pyplot.figure(figsize=(6.03,3.015))

for ivar in range(nvars):
    for i in range(ntest):

        # read data
        testdir = inputdir+"/"+testnames[i]+"/"+varnames[ivar]
        testprefix = testdir+"/"+varnames[ivar]+"_"
        test = myvars.taylor_var(testprefix,varnames[ivar])

        cntldir = inputdir+"/"+cntlnames[ivar]+"/"+varnames[ivar]
        cntlprefix = cntldir+"/"+varnames[ivar]+"_"
        cntl = myvars.taylor_var(cntlprefix,varnames[ivar])

        # regrid
        if len(test.lat) < len(cntl.lat) or len(test.lon) < len(cntl.lon):
            cntl.regrid(test.lat,test.lon)
        elif len(cntl.lat) < len(test.lat) or len(cntl.lon) < len(test.lon):
            test.regrid(cntl.lat,cntl.lon)

        # mask missing data
        test.mask_area(cntl)
        cntl.mask_area(test)

        # subset region
        test.region_subset((-60.0,60.0),(0.0,360.0))
        cntl.region_subset((-60.0,60.0),(0.0,360.0))

        # calculate taylor statistics
        ratio[i],cc[i],bias[i],rmserr[i] = test.taylor_statistics(cntl)
        
    # make figure
    ax = figure.add_subplot(1,2,ivar+1,frameon=False)
    taylor_diagram = taylor.Taylor_diagram(
            ax,ratio,cc,bias,
            casecolors=testcolors,
        )

# test labels
figure.axes[0].text(
        0.97,0.67,"CAM3",
        color="SeaGreen",
        horizontalalignment="left"
    )
figure.axes[0].text(
        0.78,0.55,"CAM4",
        color="RoyalBlue",
        horizontalalignment="right"
    )
figure.axes[0].text(
        0.93,0.50,"CAM5",
        color="Firebrick",
        horizontalalignment="left"
    )
figure.axes[0].text(
        1.21,0.26,"ERBE",
        color="0.55",
        horizontalalignment="center"
    )

# plot labels
figure.axes[0].text(
        0.05,0.15,r"Shortwave CRE",
        color="Black",
        horizontalalignment="left"
    )
figure.axes[1].text(
        0.05,0.15,r"Longwave CRE",
        color="Black",
        horizontalalignment="left"
    )
figure.axes[0].text(
        0.05,0.05,r"60S - 60N",
        color="Black",
        horizontalalignment="left",
    )
figure.axes[1].text(
        0.05,0.05,r"60S - 60N",
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
figure.savefig("../graphics/cre_camamip_taylor.pdf",format="pdf")
