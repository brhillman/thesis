import numpy as numpy
import matplotlib.pyplot as pyplot
import pylab as pylab
import myvars as myvars

# sum bins for each cloud type and annotate plot with values
def sum_cth(ax,histogram,tau,(tau1,tau2),cth,(cth1,cth2)):
    if tau1 is None:
        tind1 = 0
    else:
        tind1 = next(idx for idx,value in enumerate(tau) if value > tau1)
    if tau2 is None:
        tind2 = len(tau)
    else:
        tind2 = next(idx for idx,value in enumerate(tau) if value > tau2)
    if cth1 is None:
        cind1 = 0
    else:
        cind1 = next(idx for idx,value in enumerate(cth) if value > cth1)
    if cth2 is None:
        cind2 = len(cth)
    else:
        cind2 = next(idx for idx,value in enumerate(cth) if value > cth2)
    boxmean = histogram[cind1:cind2,tind1:tind2].sum()
    ax.text(
            (tind1+tind2)/2.0,(cind1+cind2)/2.0,"%.1f"%boxmean,
            color="black",
            horizontalalignment="center",
            verticalalignment="center"
        )
    ax.plot([tind1,tind1],[cind1,cind2],color="black",linestyle="dotted")
    ax.plot([tind2,tind2],[cind1,cind2],color="black",linestyle="dotted")
    ax.plot([tind1,tind2],[cind1,cind1],color="black",linestyle="dotted")
    ax.plot([tind1,tind2],[cind2,cind2],color="black",linestyle="dotted")


# sum bins for each cloud type and annotate plot with values
def sum_plev(ax,histogram,tau,(tau1,tau2),plev,(plev1,plev2)):
    if tau1 is None:
        tind1 = 0
    else:
        tind1 = next(idx for idx,value in enumerate(tau) if value > tau1)
    if tau2 is None:
        tind2 = len(tau)
    else:
        tind2 = next(idx for idx,value in enumerate(tau) if value > tau2)
    if plev1 is None:
        cind1 = 0
    else:
        cind1 = next(idx for idx,value in enumerate(plev) if value < plev1)
    if plev2 is None:
        cind2 = len(plev)
    else:
        cind2 = next(idx for idx,value in enumerate(plev) if value < plev2)
    boxmean = histogram[cind1:cind2,tind1:tind2].sum()
    ax.text(
            (tind1+tind2)/2.0,(cind1+cind2)/2.0,"%.1f"%boxmean,
            color="black",
            horizontalalignment="center",
            verticalalignment="center"
        )
    ax.plot([tind1,tind1],[cind1,cind2],color="black",linestyle="dotted")
    ax.plot([tind2,tind2],[cind1,cind2],color="black",linestyle="dotted")
    ax.plot([tind1,tind2],[cind1,cind1],color="black",linestyle="dotted")
    ax.plot([tind1,tind2],[cind2,cind2],color="black",linestyle="dotted")


# fix fonts
params = {'axes.titlesize': 10,
          'axes.labelsize': 10,
          'text.fontsize': 10,
          'legend.fontsize': 10,
          'xtick.labelsize': 10,
          'ytick.labelsize': 10,
          'font.family': 'serif',
          'font.family': 'serif',
          'font.size': 10,
          'text.usetex': True}
pylab.rcParams.update(params)

basedir = "~/cmor/climo"
maskfile = basedir+"/misr/clMISR/clMISR_ANN.nc"

# test and control variables
testinfo = numpy.array([
        ["clMISR" ,"cam3_hot4k","CAM3 HOT4K"],
        ["clMISR" ,"am2_hot4k" ,"AM2 HOT4K" ],
        ["clisccp","cam3_hot4k","CAM3 HOT4K"],
        ["clisccp","am2_hot4k" ,"AM2 HOT4K" ],
        ["clmodis","cam3_hot4k","CAM3 HOT4K"],
        ["clmodis","am2_hot4k" ,"AM2 HOT4K" ],
    ])
cntlinfo = numpy.array([
        ["clMISR" ,"cam3_amip","CAM3 AMIP"],
        ["clMISR" ,"am2_amip" ,"AM2 AMIP" ],
        ["clisccp","cam3_amip","CAM3 AMIP"],
        ["clisccp","am2_amip" ,"AM2 AMIP" ],
        ["clmodis","cam3_amip","CAM3 AMIP"],
        ["clmodis","am2_amip" ,"AM2 AMIP" ],
    ])
testvars = testinfo[:,0]
testnames = testinfo[:,1]
testlabels = testinfo[:,2]
cntlvars = cntlinfo[:,0]
cntlnames = cntlinfo[:,1]
cntllabels = cntlinfo[:,2]
nvars = len(testvars)

# region boundaries
lat1 = 15.0
lat2 = 35.0
lon1 = 220.0
lon2 = 250.0

# open figure
figure = pyplot.figure(figsize=(6.03,6.03))

# contour levels
cmap = pyplot.cm.get_cmap("RdBu_r",20)
clevels = numpy.arange(0,5.1,20)

# axes bounds
tau_bnds = numpy.array([0.3,1.3,3.6,9.4,23.0,60.0,379.0])
plev_bnds = numpy.array([1100.0,800.0,680.0,560.0,440.0,310.0,180.0,0.0])
cth_bnds = numpy.array([
        0.0,0.5,1.0,1.5,2.0,2.5,3.0,4.0,5.0,7.0,9.0,11.0,13.0,15.0,17.0,23.0
    ])

for i in range(nvars):

    # read joint histogram
    inputdir = basedir+"/"+testnames[i]+"/"+testvars[i]
    filename = inputdir+"/"+testvars[i]+"_JJA.nc"
    test = myvars.hist_var(filename,testvars[i])

    inputdir = basedir+"/"+cntlnames[i]+"/"+cntlvars[i]
    filename = inputdir+"/"+cntlvars[i]+"_JJA.nc"
    cntl = myvars.hist_var(filename,cntlvars[i])

    # need to mask land here for misr
    if testvars[i] == "clMISR":
        maskvar = myvars.hist_var(maskfile,"clMISR")
        maskvar.regrid(test.lat,test.lon)
        test.mask_area(maskvar)

    # average over region
    test.region_subset((lat1,lat2),(lon1,lon2))
    cntl.region_subset((lat1,lat2),(lon1,lon2))

    # calculate histogram
    histogram = test.weighted_area_average() - cntl.weighted_area_average()

    # clMISR has tau,cth, we want cth,tau
    if testvars[i] == "clMISR":
        histogram = numpy.transpose(histogram)

    # remove bins below minimum cloud optical thickness
    tind = next(idx for idx,tau in enumerate(test.tau) if tau > tau_bnds[0])
    histogram = histogram[:,tind:]
    tau = test.tau[tind:]

    # remove bins below minimum cloud top height
    if testvars[i] == "clMISR":
        hind = next(idx for idx,cth in enumerate(test.cth) if cth > cth_bnds[0])
        histogram = histogram[hind:,:]
        cth = test.cth[hind:]
    else:
        plev = test.plev

    # plot
    ax = figure.add_subplot(3,2,i+1)
    plot = ax.pcolor(
            histogram,
            cmap=cmap,
            vmin=-3.0,vmax=3.0,
        )

    # axis tick labels off by default
    ax.set_xticklabels([])
    ax.set_yticklabels([])

    # plot title
    mean = histogram.sum()
    if testvars[i] == "clMISR":
        ax.set_title(testlabels[i]+" - "+cntllabels[i]+" (MISR)")
    elif testvars[i] == "clisccp":
        ax.set_title(testlabels[i]+" - "+cntllabels[i]+" (ISCCP)")
    elif testvars[i] == "clmodis":
        ax.set_title(testlabels[i]+" - "+cntllabels[i]+" (MODIS)")

    # sum bins
    if testvars[i] == "clMISR":
        sum_cth(ax,histogram,tau,(0.3,3.6),cth,(0.0,3000.0))
        sum_cth(ax,histogram,tau,(0.3,3.6),cth,(3000.0,7000.0))
        sum_cth(ax,histogram,tau,(0.3,3.6),cth,(7000.0,None))
        sum_cth(ax,histogram,tau,(3.6,23.0),cth,(0.0,3000.0))
        sum_cth(ax,histogram,tau,(3.6,23.0),cth,(3000.0,7000.0))
        sum_cth(ax,histogram,tau,(3.6,23.0),cth,(7000.0,None))
        sum_cth(ax,histogram,tau,(23.0,None),cth,(0.0,3000.0))
        sum_cth(ax,histogram,tau,(23.0,None),cth,(3000.0,7000.0))
        sum_cth(ax,histogram,tau,(23.0,None),cth,(7000.0,None))
    else:
        sum_plev(ax,histogram,tau,(0.3,3.6),plev,(None,68000.0))
        sum_plev(ax,histogram,tau,(0.3,3.6),plev,(68000.0,44000.0))
        sum_plev(ax,histogram,tau,(0.3,3.6),plev,(44000.0,None))
        sum_plev(ax,histogram,tau,(3.6,23.0),plev,(None,68000.0))
        sum_plev(ax,histogram,tau,(3.6,23.0),plev,(68000.0,44000.0))
        sum_plev(ax,histogram,tau,(3.6,23.0),plev,(44000.0,None))
        sum_plev(ax,histogram,tau,(23.0,None),plev,(None,68000.0))
        sum_plev(ax,histogram,tau,(23.0,None),plev,(68000.0,44000.0))
        sum_plev(ax,histogram,tau,(23.0,None),plev,(44000.0,None))
    ax.set_xlim([0,len(tau)])

    # set axis limits
    if testvars[i] == "clMISR":
        ax.set_ylim([0,len(cth)])
    else:
        ax.set_ylim([0,len(plev)])

# optical thickness labels
xticklabels = []
for i in range(len(tau_bnds)):
    if tau_bnds[i] < 10.0:
        xticklabels.append("%.1f"%tau_bnds[i])
    else:
        xticklabels.append("%.0f"%tau_bnds[i])
for i in [4,5]:
    ax = figure.add_axes(figure.axes[i])
    ax.set_xticklabels(xticklabels)
    ax.set_xlabel("Cloud optical thickness")

# cloud top height labels for the first row
cthticklabels = []
cthticks = numpy.arange(len(cth_bnds))
for i in range(len(cth_bnds)):
    cthticklabels.append("%.1f"%cth_bnds[i])
ax = figure.add_axes(figure.axes[0])
ax.set_yticks(cthticks[::2])
ax.set_yticklabels(cthticklabels[::2])
ax.set_ylabel("CTH (km)")
ax.yaxis.set_label_coords(-0.2,0.5)

# cloud top pressure labels for the second and third rows
plevticklabels = []
for i in range(len(plev_bnds)):
    plevticklabels.append("%.0f"%plev_bnds[i])
for i in [2,4]:
    ax = figure.add_axes(figure.axes[i])
    ax.set_yticklabels(plevticklabels)
    #ax.set_ylabel("Cloud top pressure (hPa)")
    ax.set_ylabel("CTP (hPa)")
    ax.yaxis.set_label_coords(-0.2,0.5)

# adjust subplots
figure.subplots_adjust(left   = 0.10)
figure.subplots_adjust(bottom = 0.20)
figure.subplots_adjust(top    = 0.95)
figure.subplots_adjust(right  = 0.95)
figure.subplots_adjust(hspace = 0.20)
figure.subplots_adjust(wspace = 0.20)

# common colorbar
cax = figure.add_axes([0.10,0.08,0.85,0.03])
cb = pyplot.colorbar(
        plot,
        orientation='horizontal',
        cax=cax,
        drawedges=False,
    )
cb.set_label(r"Frequency (\%)")

# save figure
figure.savefig("../graphics/hist2d_cmip3hot4k_california.pdf",format="pdf")
