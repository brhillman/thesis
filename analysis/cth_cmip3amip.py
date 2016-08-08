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

basedir = "~/cmor/climo"

# test variables
testinfo = numpy.array([
        ["cam3_amip","CAM3","SeaGreen"],
        ["am2_amip","AM2","Coral"]
    ])
testnames  = testinfo[:,0]
testlabels = testinfo[:,1]
testcolors = testinfo[:,2]
ntest = len(testnames)

# control variable
varname = "clMISR"
cntlname = "misr"
cntllabel = "MISR"

# open figure
golden = (1+numpy.sqrt(5.0))/2.0
width = 6.03
height = width/golden
figure = pyplot.figure(figsize=(width,height))
ax = figure.add_subplot(111)

# read control data
cntlfile = basedir+"/"+cntlname+"/"+varname+"/"+varname+"_ANN.nc"
cntl = myvars.cth_var(cntlfile,varname)

# plot control data
cntl.calculate()
ax.plot(
        cntl.histogram,range(len(cntl.histogram)),'.',
        color="black",
        linestyle="solid",
        linewidth=2.0,
        markersize=10.0,
        label=cntllabel,
    )

for i in range(ntest):

    # read test data
    testfile = basedir+"/"+testnames[i]+"/"+varname+"/"+varname+"_ANN.nc"
    test = myvars.cth_var(testfile,varname)

    # mask ice and land
    maskvar = myvars.cth_var(cntlfile,varname)
    maskvar.regrid(test.lat,test.lon)
    test.mask_area(maskvar)

    # plot test data
    test.calculate()
    ax.plot(
            test.histogram,range(len(test.histogram)),'.',
            color=testcolors[i],
            linestyle="solid",
            linewidth=2.0,
            markersize=10.0,
            label=testlabels[i]
        )

yticklabels = []
for i in range(test.cth_bnds.shape[0]):
    yticklabels.append(test.cth_bnds[i,0]/1000.0)

ax.set_yticks(numpy.array(range(len(yticklabels)))[::2]-0.5)
ax.set_yticklabels(yticklabels[::2])
ax.set_xlabel("Frequency (\%)")
ax.set_ylabel("Cloud top height (km)")
ax.grid()

# add a legend
pyplot.legend(loc="upper right",handlelength=3)

# adjust
figure.subplots_adjust(left=0.125)
figure.subplots_adjust(bottom=0.15)
figure.subplots_adjust(top=0.95)
figure.subplots_adjust(right=0.95)

# save figure
figure.savefig("../graphics/cth_cmip3amip.pdf",format="pdf")
