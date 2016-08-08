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
maskfile = basedir+"/misr/clMISR/clMISR_ANN.nc"

# variable info
varinfo = numpy.array([
        ["clMISR" ,"misr"     ,"MISR" ],
        ["clMISR" ,"cam3_amip","CAM3" ],
        ["clMISR" ,"am2_amip" ,"AM2"  ],
        ["clisccp","isccp"    ,"ISCCP"],
        ["clisccp","cam3_amip","CAM3" ],
        ["clisccp","am2_amip" ,"AM2"  ],
        ["clmodis","modis"    ,"MODIS"],
        ["clmodis","cam3_amip","CAM3" ],
        ["clmodis","am2_amip" ,"AM2"  ],
    ])
varnames = varinfo[:,0]
casenames = varinfo[:,1]
caselabels = varinfo[:,2]
nvars = len(varnames)

# region boundaries
lat1 = 15.0
lat2 = 35.0
lon1 = 150.0
lon2 = 210.0

# open figure
figure = pyplot.figure(figsize=(6.03,6.03))

# contour levels
cmap = pyplot.cm.get_cmap("Blues",20)
clevels = numpy.arange(0,5.1,20)

# axes bounds
tau_bnds = numpy.array([0.3,1.3,3.6,9.4,23.0,60.0,379.0])
plev_bnds = numpy.array([1100.0,800.0,680.0,560.0,440.0,310.0,180.0,0.0])
cth_bnds = numpy.array([
        0.0,0.5,1.0,1.5,2.0,2.5,3.0,4.0,5.0,7.0,9.0,11.0,13.0,15.0,17.0,23.0
    ])

for i in range(nvars):

    # read joint histogram
    inputdir = basedir+"/"+casenames[i]+"/"+varnames[i]
    filename = inputdir+"/"+varnames[i]+"_JJA.nc"
    test = myvars.hist_var(filename,varnames[i])

    # need to mask land here for misr
    if varnames[i] == "clMISR":
        maskvar = myvars.hist_var(maskfile,"clMISR")
        maskvar.regrid(test.lat,test.lon)
        test.mask_area(maskvar)

    # average over region
    test.region_subset((lat1,lat2),(lon1,lon2))
    histogram = test.weighted_area_average()

    # clMISR has tau,cth, we want cth,tau
    if varnames[i] == "clMISR":
        histogram = numpy.transpose(histogram)

    # remove bins below minimum cloud optical thickness
    tind = next(idx for idx,tau in enumerate(test.tau) if tau > tau_bnds[0])
    histogram = histogram[:,tind:]
    tau = test.tau[tind:]

    # remove bins below minimum cloud top height
    if varnames[i] == "clMISR":
        hind = next(idx for idx,cth in enumerate(test.cth) if cth > cth_bnds[0])
        histogram = histogram[hind:,:]
        cth = test.cth[hind:]
    else:
        plev = test.plev

    # plot
    ax = figure.add_subplot(3,3,i+1)
    plot = ax.pcolor(
            histogram,
            cmap=cmap,
            vmin=0.0,vmax=5.0,
        )

    # axis tick labels off by default
    ax.set_xticklabels([])
    ax.set_yticklabels([])

    # plot title
    mean = histogram.sum()
    if varnames[i] == "clMISR":
        if caselabels[i] != "MISR":
            ax.set_title(caselabels[i]+r" MISR SIM (%.0f%s)"%(mean,"\%"))
        else:
            ax.set_title(caselabels[i]+r" (%.0f%s)"%(mean,"\%"))
    elif varnames[i] == "clisccp":
        if caselabels[i] != "ISCCP":
            ax.set_title(caselabels[i]+r" ISCCP SIM (%.0f%s)"%(mean,"\%"))
        else:
            ax.set_title(caselabels[i]+r" (%.0f%s)"%(mean,"\%"))
    elif varnames[i] == "clmodis":
        if caselabels[i] != "MODIS":
            ax.set_title(caselabels[i]+r" MODIS SIM (%.0f%s)"%(mean,"\%"))
        else:
            ax.set_title(caselabels[i]+r" (%.0f%s)"%(mean,"\%"))

    # sum bins
    if varnames[i] == "clMISR":
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
    if varnames[i] == "clMISR":
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
for i in [6,7,8]:
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
ax.yaxis.set_label_coords(-0.3,0.5)

# cloud top pressure labels for the second and third rows
plevticklabels = []
for i in range(len(plev_bnds)):
    plevticklabels.append("%.0f"%plev_bnds[i])
for i in [3,6]:
    ax = figure.add_axes(figure.axes[i])
    ax.set_yticklabels(plevticklabels)
    #ax.set_ylabel("Cloud top pressure (hPa)")
    ax.set_ylabel("CTP (hPa)")
    ax.yaxis.set_label_coords(-0.3,0.5)

# adjust subplots
figure.subplots_adjust(left   = 0.15)
figure.subplots_adjust(bottom = 0.20)
figure.subplots_adjust(top    = 0.95)
figure.subplots_adjust(right  = 0.95)
figure.subplots_adjust(hspace = 0.20)
figure.subplots_adjust(wspace = 0.20)

# common colorbar
cax = figure.add_axes([0.15,0.08,0.80,0.03])
cb = pyplot.colorbar(
        plot,
        orientation='horizontal',
        cax=cax,
        drawedges=False,
    )
cb.set_label(r"Frequency (\%)")

# save figure
figure.savefig("../graphics/hist2d_cmip3amip_hawaiian.pdf",format="pdf")
