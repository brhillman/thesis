import numpy as numpy
import matplotlib as matplotlib
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

basedir = "~/cmor/climo"
maskfile = basedir+"/misr/cltmisr/cltmisr_ANN.nc"

varinfo = numpy.array([
        [
                ["cltisccp","ISCCP ","total"],
                ["cltmisr","MISR ","total"],
                ["cltmodis","MODIS ","total"],
            ],
        [
                ["cllisccp","ISCCP","low"],
                ["cllmisr","MISR","low"],
                ["cllmodis","MODIS","low"],
            ],
        [
                ["clmisccp","ISCCP","mid"],
                ["clmmisr","MISR","mid"],
                ["clmmodis","MODIS","mid"],
            ],
        [
                ["clhisccp","ISCCP","high"],
                ["clhmisr","MISR","high"],
                ["clhmodis","MODIS","high"],
            ]
    ])

# we will want to mask these variables
misrvars = ["cltmisr","cllmisr","clmmisr","clhmisr"]

# we want to use the cloud mask variables for MODIS
modisvars = ["cltmodis","cllmodis","clmmodis","clhmodis"]

caseinfo = numpy.array([
        ["obs","retrievals","0.55"],
        ["cam3_amip","CAM3","SeaGreen"],
        ["cam4_1deg_release_amip","CAM4","RoyalBlue"],
        ["cam5_1deg_release_amip","CAM5","Firebrick"],
    ])
varnames = varinfo[:,:,0]
varlabels = varinfo[:,:,1]
vartitles = varinfo[:,:,2]

casenames = caseinfo[:,0]
caselabels = caseinfo[:,1]
casecolors = caseinfo[:,2]

nplots = len(varnames[:,0])
nvars = len(varnames[0,:])
ncases = len(casenames)

width = 6.03
height = 4.5
figure = pyplot.figure(figsize=(width,height))

for iplot in range(nplots):
    ax = figure.add_subplot(2,2,iplot+1)
    for ivar in range(nvars):
        for icase in range(ncases):

            # for MODIS we use the cloud mask variables
            if varnames[iplot,ivar] in modisvars:
                if casenames[icase] == "obs":
                    varname = varnames[iplot,ivar]+"_mask"
                else:
                    varname = varnames[iplot,ivar]
            else:
                varname = varnames[iplot,ivar]

            # read data
            testdir = basedir+"/"+casenames[icase]+"/"+varname
            testfile = testdir+"/"+varname+"_ANN.nc"
            test = myvars.map_var(filename=testfile,varname=varname)

            # mask ice and land for MISR variables
            if varname in misrvars:
                maskvar = myvars.map_var(filename=maskfile,varname="cltmisr")
                maskvar.regrid(test.lat,test.lon)
                test.mask_area(maskvar)

            # calculate mean and plot
            mean = test.weighted_area_average()
            if iplot == 1 and ivar == 0:
                ax.bar(icase+ivar*(ncases+1)-ncases/2.0,mean,
                        width=1.0,
                        color=casecolors[icase],
                        alpha=0.75,
                        label=caselabels[icase]
                    )
            else:
                ax.bar(icase+ivar*(ncases+1)-ncases/2.0,mean,
                        width=1.0,
                        color=casecolors[icase],
                        alpha=0.75,
                    )

            # add a legend
            if iplot == 1:
                ax.set_ylim([0.0,80.0])
                pyplot.legend(loc="upper right",handlelength=2)

    ax.set_xticks(range(ncases*(nvars))[::ncases+1])
    ax.set_xticklabels([])

# add ylabels
for idx in [0,2]:
    ax = figure.axes[idx]
    ax.set_ylabel(r"Cloud area fraction (\%)")

# add xlabels
for idx in [2,3]:
    ax = figure.axes[idx]
    ax.set_xticklabels(varlabels[0,:])

# add titles
ax = figure.axes[0]
ax.set_title("Total cloud")

ax = figure.axes[1]
ax.set_title("Low-topped cloud")

ax = figure.axes[2]
ax.set_title("Mid-topped cloud")

ax = figure.axes[3]
ax.set_title("High-topped cloud")

# adjust subplots
figure.subplots_adjust(left=0.08)
figure.subplots_adjust(bottom=0.05)
figure.subplots_adjust(top=0.95)
figure.subplots_adjust(right=0.98)
figure.subplots_adjust(wspace=0.20)
figure.subplots_adjust(hspace=0.20)

# save figure
figure.savefig("../graphics/cldtypes_camamip_bar.pdf",format="pdf")
