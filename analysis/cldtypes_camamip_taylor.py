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

def make_plot(figure,varnames,cntlnames,cntllabels,testnames,testcolors,iplot):
    # initialize statistics variables
    nvars = len(varnames)
    ntest = len(testnames)
    ratio  = numpy.zeros([nvars,ntest])
    cc     = numpy.zeros([nvars,ntest])
    bias   = numpy.zeros([nvars,ntest])
    rmserr = numpy.zeros([nvars,ntest])

    # loop over variables and test cases
    for ivar in range(nvars):
        for itest in range(ntest):

            # test variable name
            testvarname = varnames[ivar]

            # for MODIS variables we will use the cloud mask
            if cntlnames[ivar] == "modis":
                cntlvarname = varnames[ivar]+"_mask"
            else:
                cntlvarname = varnames[ivar]

            # read data
            testdir = inputdir+"/"+testnames[itest]+"/"+testvarname
            testprefix = testdir+"/"+testvarname+"_"
            test = myvars.taylor_var(testprefix,testvarname)

            cntldir = inputdir+"/"+cntlnames[ivar]+"/"+cntlvarname
            cntlprefix = cntldir+"/"+cntlvarname+"_"
            cntl = myvars.taylor_var(cntlprefix,cntlvarname)

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
            ratio0,cc0,bias0,rmserr0 = test.taylor_statistics(cntl)
            ratio[ivar,itest] = ratio0
            cc[ivar,itest] = cc0
            bias[ivar,itest] = bias0
            rmserr[ivar,itest] = rmserr0

    # make figure
    ax = figure.add_subplot(2,2,iplot+1,frameon=False)
    taylor_diagram = taylor.Taylor_diagram(
            ax,ratio,cc,bias,
            casecolors=testcolors,
            varlabels=cntllabels
        )


inputdir = "~/cmor/climo"
testinfo = numpy.array([
        ["cam3_amip","CAM3","SeaGreen"],
        ["cam4_1deg_release_amip","CAM4","RoyalBlue"],
        ["cam5_1deg_release_amip","CAM5","Firebrick"],
    ])
testnames = testinfo[:,0]
testlabels = testinfo[:,1]
testcolors = testinfo[:,2]

cntlnames = ["isccp","misr","modis"]
cntllabels = ["A","B","C"]

# open figure
figure = pyplot.figure(figsize=(6.03,6.50))

# make plots
varnames = ["cltisccp","cltmisr","cltmodis"]
make_plot(figure,varnames,cntlnames,cntllabels,testnames,testcolors,0)

varnames = ["cllisccp","cllmisr","cllmodis"]
make_plot(figure,varnames,cntlnames,cntllabels,testnames,testcolors,1)

varnames = ["clmisccp","clmmisr","clmmodis"]
make_plot(figure,varnames,cntlnames,cntllabels,testnames,testcolors,2)

varnames = ["clhisccp","clhmisr","clhmodis"]
make_plot(figure,varnames,cntlnames,cntllabels,testnames,testcolors,3)

# test labels
figure.axes[0].text(
        0.55,1.02,"CAM3",
        color="SeaGreen",
        horizontalalignment="right"
    )
figure.axes[0].text(
        0.42,0.63,"CAM4",
        color="RoyalBlue",
        horizontalalignment="right"
    )
figure.axes[0].text(
        1.00,0.45,"CAM5",
        color="Firebrick",
        horizontalalignment="left"
    )

# control labels
figure.axes[0].text(
        0.05,1.35,r"A - ISCCP",
        color="Black",
        horizontalalignment="left"
    )
figure.axes[0].text(
        0.05,1.25,r"B - MISR",
        color="Black",
        horizontalalignment="left"
    )
figure.axes[0].text(
        0.05,1.15,r"C - MODIS",
        color="Black",
        horizontalalignment="left"
    )

# plot labels
figure.axes[0].text(
        0.05,0.15,r"Total cloud",
        color="Black",
        horizontalalignment="left"
    )
figure.axes[1].text(
        0.05,0.15,r"Low-topped cloud",
        color="Black",
        horizontalalignment="left"
    )
figure.axes[2].text(
        0.05,0.15,r"Mid-topped cloud",
        color="Black",
        horizontalalignment="left"
    )
figure.axes[3].text(
        0.05,0.15,r"High-topped cloud",
        color="Black",
        horizontalalignment="left"
    )
for i in range(4):
    figure.axes[i].text(
            0.05,0.05,r"60S - 60N",
            color="Black",
            horizontalalignment="left",
        )

# adjust subplots
figure.subplots_adjust(left   = 0.01)
figure.subplots_adjust(bottom = 0.05)
figure.subplots_adjust(top    = 1.00)
figure.subplots_adjust(right  = 0.95)
figure.subplots_adjust(hspace = 0.05)
figure.subplots_adjust(wspace = 0.10)
 
# save figure
figure.savefig("../graphics/cldtypes_camamip_taylor.pdf",format="pdf")
