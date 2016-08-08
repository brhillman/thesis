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
maskfile = inputdir+"/misr/clMISR/clMISR_ANN.nc"
  
varinfo  = numpy.array([
        ["clisccp","solid" ,"isccp","ISCCP","black"],
        ["clmodis","dashed","modis","MODIS","black"],
        ["clMISR","dotted","misr","MISR" ,"black"]
    ])
testinfo = numpy.array([
        ["cam3_amip","CAM3 T42","SeaGreen"],
        ["cam3_amip_fv","CAM3 FV","MediumAquamarine"],
        ["cam4_1deg_release_amip","CAM4","RoyalBlue"],
        ["cam5_1deg_release_amip","CAM5","Firebrick"],
    ])
varnames = varinfo[:,0]
varlines = varinfo[:,1]
cntlnames = varinfo[:,2]
cntllabels = varinfo[:,3]
cntlcolors = varinfo[:,4]

testnames = testinfo[:,0]
testlabels = testinfo[:,1]
testcolors = testinfo[:,2]
ntest = len(testnames)
nvars = len(varnames)

# axes bounds
tau_bnds = [1.3,3.6,9.4,23.0,60.0,379.0]
histogram = numpy.zeros(len(tau_bnds))

# open figure
golden = (1+numpy.sqrt(5.0))/2.0
width = 6.03
height = width/golden
figure = pyplot.figure(figsize=(width,height))
ax = figure.add_subplot(111)

# we will want to save each curve so we can plot between the min and max
testdata = numpy.zeros([nvars,ntest,len(tau_bnds)-1])
cntldata = numpy.zeros([nvars,len(tau_bnds)-1])

for ivar in range(nvars):

    # read control data
    directory = inputdir+"/"+cntlnames[ivar]+"/"+varnames[ivar]
    filename = directory+"/"+varnames[ivar]+"_ANN.nc"
    cntl = myvars.tau_var(filename,varnames[ivar])

    # apply additional mask
    maskvar = myvars.tau_var(maskfile,"clMISR")
    maskvar.regrid(cntl.lat,cntl.lon)
    cntl.mask_area(maskvar)

    # sum histogram
    cntl.calculate()
    cntldata[ivar,:] = cntl.histogram

    # plot control data
    ax.plot(
            cntl.histogram,'.',
            color=cntlcolors[ivar],
            linestyle=varlines[ivar],
            linewidth=2.0,
            markersize=10.0,
            label=cntllabels[ivar]
        )

    for itest in range(ntest):

        # read test data
        directory = inputdir+"/"+testnames[itest]+"/"+varnames[ivar]
        filename = directory+"/"+varnames[ivar]+"_ANN.nc"
        test = myvars.tau_var(filename,varnames[ivar])

        # apply additional mask
        maskvar = myvars.tau_var(maskfile,"clMISR")
        maskvar.regrid(test.lat,test.lon)
        test.mask_area(maskvar)

        # sum histogram
        test.calculate()
        testdata[ivar,itest,:] = test.histogram

        # plot test data
        ax.plot(
                test.histogram,'.',
                color=testcolors[itest],
                linestyle=varlines[ivar],
                linewidth=2.0,
                markersize=10.0,
            )

# shade between curves
mintau = numpy.zeros(len(tau_bnds)-1)
maxtau = numpy.zeros(len(tau_bnds)-1)
for i in range(len(tau_bnds)-1):
    mintau[i] = numpy.min(cntldata[:,i])
    maxtau[i] = numpy.max(cntldata[:,i])
ax.fill_between(
        range(len(tau_bnds)-1),mintau,maxtau,
        facecolor="black",
        edgecolor="none",
        alpha=0.25
    )
for itest in range(ntest):
    for i in range(len(tau_bnds)-1):
        mintau[i] = numpy.min(testdata[:,itest,i])
        maxtau[i] = numpy.max(testdata[:,itest,i])
    for ivar in range(nvars-1):
        ax.fill_between(
                range(len(tau_bnds)-1),mintau,maxtau,
                facecolor=testcolors[itest],
                edgecolor="none",
                alpha=0.25
            )

# set axes ticks and labels
ax.set_xticks(numpy.array(range(len(tau_bnds)))-0.5)
ax.set_xticklabels(tau_bnds)
ax.set_xlabel("Cloud optical thickness")
ax.set_ylabel("Frequency (\%)")
ax.grid()

# add a legend
pyplot.legend(loc="lower center",handlelength=3)

# add some annotations
ax.text(0.45,19.0,"retrievals",color="black",horizontalalignment="left")
ax.text(3.40,13.0,"CAM3",color="SeaGreen",horizontalalignment="left") 
ax.text(3.10,15.0,"CAM3 FV",color="MediumAquamarine",horizontalalignment="left")
ax.text(3.70, 5.8,"CAM4",color="RoyalBlue",horizontalalignment="left")
ax.text(2.10,17.1,"CAM5",color="Firebrick",horizontalalignment="left")

# adjust subplots
figure.subplots_adjust(left=0.1)
figure.subplots_adjust(bottom=0.15)
figure.subplots_adjust(top=0.95)
figure.subplots_adjust(right=0.95)

# save figure
figure.savefig("../graphics/tau_camamip.pdf",format="pdf")
