import numpy as numpy
import Scientific.IO.NetCDF as netcdf
import matplotlib as matplotlib
import matplotlib.pyplot as pyplot
from mpl_toolkits.basemap import Basemap
import pylab as pylab
import myvars

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

input_dir = "~/cmor/climo"

caseinfo = numpy.array([
        ["cltisccp"     ,"isccp","ISCCP"           ],
        ["cltmisr"      ,"misr" ,"MISR"            ],
        ["cltmodis"     ,"modis","MODIS retrievals"],
        ["cltmodis_mask","modis","MODIS cloud mask"]
    ])
varnames = caseinfo[:,0]
casenames = caseinfo[:,1]
caselabels = caseinfo[:,2]
ncases = len(varnames)

figure = pyplot.figure(figsize=(6.03,4.35))
for i in range(ncases):

    # read data
    test_dir = input_dir+"/"+casenames[i]+"/"+varnames[i]
    testfile = test_dir+"/"+varnames[i]+"_ANN.nc"
    test = myvars.map_var(filename=testfile,varname=varnames[i])

    # calculate means and standard deviation
    mean = test.weighted_area_average()

    # make figure
    ax = figure.add_subplot(2,2,i+1)

    # transform coordinates
    datamap = Basemap(projection='robin',lon_0=0.5*(test.lon[0]+test.lon[-1]))
    x,y = datamap(*numpy.meshgrid(test.lon,test.lat))

    # plot data
    cmap = pyplot.cm.Blues
    clevels = numpy.arange(5,96,5)
    plot  = datamap.contourf(x,y,test.data,clevels,
                             extend="both",cmap=cmap)
    pyplot.title(caselabels[i]+" (%.0f%s)"%(mean,"\%"))

    # draw lines
    datamap.drawcoastlines(linewidth=0.5)
    datamap.drawmapboundary(linewidth=0.5)

# adjust subplots
figure.subplots_adjust(left   = 0.00)
figure.subplots_adjust(bottom = 0.12)
figure.subplots_adjust(top    = 0.99)
figure.subplots_adjust(right  = 1.00)
figure.subplots_adjust(hspace = 0.01)
figure.subplots_adjust(wspace = 0.01)

# common colorbar
ticks = [-40,-20,0,20,40]
cax = figure.add_axes([0.10,0.10,0.8,0.04])
cb = pyplot.colorbar(plot,
                     orientation='horizontal',cax=cax,
                     drawedges=False)
cb.set_label(r"Cloud area fraction (\%)")

# save figure
figure.savefig("../graphics/clt_retrievals_map.pdf",format="pdf")
