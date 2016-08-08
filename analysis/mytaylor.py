import numpy
import matplotlib.artist as artist
import matplotlib.patches as patches
import matplotlib.lines as lines
import matplotlib.pylab as pylab
from matplotlib.axes import Axes
import Scientific.IO.NetCDF as netcdf

class Taylor_diagram():
    """
    Custom axes class to facilitate making Taylor diagrams. Adds the taylor
    method to draw the diagram.
    """

    def __init__(self,ax,ratio,cc,bias,rms=None,casecolors=None,varlabels=None):

        # default casecolors
        if casecolors is None:
            casecolors = [
                    'blue','red','green','cyan','magenta','yellow',
                    '0.75','0.50','0.25','0.10'
                ]

        # copy axes instance
        self.ax = ax

        # plot size
        self.xymax = 1.50

        # draw axes
        self.draw_axes()

        # add some reference arcs
        self.draw_ratio(1.0,linestyle='solid',color='black')
        self.draw_ratio(numpy.min(ratio),linestyle='dotted',color='black')
        self.draw_ratio(numpy.max(ratio),linestyle='dotted',color='black')

        # draw some reference lines
        self.draw_radii(numpy.min(cc))
        self.draw_radii(numpy.max(cc))

        # draw points
        self.draw_point(
                ratio,cc,bias,
                bubblecolors=casecolors,
                varlabels=varlabels,
            )

        # draw rms circles
        if rms is not None:
            self.draw_rms(min(rms))
            self.draw_rms(max(rms))

    def polar_transform(self,r,theta):
        x = r*numpy.cos(theta)
        y = r*numpy.sin(theta)
        return x,y

    def draw_radii(self,cc,color='black',linestyle='dotted'):
        theta = numpy.arccos(cc)
        x,y = self.polar_transform(self.xymax,theta)
        self.ax.plot([0,x],[0,y],':k',linewidth=0.5)

    def draw_ratio(self,r,color='black',linestyle='dotted'):
        arc = patches.Arc(
                [0,0],2*r,2*r,
                theta1=0.0,
                theta2=90.0,
                color=color,
                linestyle=linestyle,
                linewidth=0.5,
            )
        self.ax.add_patch(arc)

    def draw_axes(self):

        # draw axes
        wedge = patches.Wedge(
                [0,0],self.xymax,
                theta1=0.0,
                theta2=90.0,
                color='black',
                fill=False
            )
        self.ax.add_patch(wedge)

        # correlation axis label
        self.ax.text(
                self.xymax*0.8,self.xymax*0.8,'Correlation',
                color='black',
                rotation='-45',
                horizontalalignment='center',
                verticalalignment='center'
            )

        # x-axis label
        self.ax.xaxis.set_label_text('Relative standard deviation')

        # x-axis tickmarks
        xmajorticks = [0.00,0.25,0.50,0.75,1.0,1.25,1.50]
        xmajorticklabels = ["0","","0.5","","1.0","","1.5"]

        # y-axis tickmarks
        ymajorticks = xmajorticks
        ymajorticklabels = []

        # correlation tickmarks
        ccmajorticks = [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9]
        ccminorticks = [0.91,0.92,0.93,0.94,0.95,0.96,0.97,0.98,0.99]
        ccmajorticklabels = [
            "","","","","0.5","","0.7","","0.9"
        ]

        # set and draw tickmarks
        self.ax.set_aspect('equal')
        self.ax.xaxis.set_ticks_position('bottom')
        self.ax.yaxis.set_ticks_position('left')

        self.ax.xaxis.set_ticks(xmajorticks)
        self.ax.xaxis.set_ticklabels(xmajorticklabels)
        self.ax.yaxis.set_ticks(ymajorticks)
        self.ax.yaxis.set_ticklabels(ymajorticklabels)


        # draw correlation tickmarks
        majorticklength = self.xymax*0.02
        minorticklength = self.xymax*0.01
        for cc in ccmajorticks:
            r = [self.xymax-majorticklength,self.xymax            ]
            theta = [numpy.arccos(cc),numpy.arccos(cc)]
            x,y = self.polar_transform(r,theta)
            self.ax.plot(x,y,'-k')
        for i,cc in enumerate(ccmajorticklabels):
            x,y = self.polar_transform(self.xymax,numpy.arccos(ccmajorticks[i]))
            self.ax.text(x,y,cc,color="black")
        for cc in ccminorticks:
            r = [self.xymax-minorticklength,self.xymax            ]
            theta = [numpy.arccos(cc),numpy.arccos(cc)]
            x,y = self.polar_transform(r,theta)
            self.ax.plot(x,y,'-k')
 
        # re-draw xy ticks to be consistent with cc ticks (not elegant)
        for xx in xmajorticks:
            x = [xx,xx                         ]
            y = [0 ,majorticklength]
            self.ax.plot(x,y,'-k')
        for yy in ymajorticks:
            x = [0 ,majorticklength]
            y = [yy,yy                         ]
            self.ax.plot(x,y,'-k')


    def draw_point(self,ratio,cc,bias,bubblecolors=None,varlabels=None,labeloffset=0.025):
        
        if len(ratio.shape) == 2:
            nvars = len(ratio[:,0])
            ncases = len(ratio[0,:])
            for ivar in range(nvars):

                # transform coordinates
                r = ratio[ivar,:]
                theta = numpy.arccos(cc[ivar,:])
                size = numpy.abs(bias[ivar,:])/2.0
                x,y = self.polar_transform(r,theta)
                for icase in range(ncases):

                    # draw "bubbles" around points
                    if bubblecolors is None: 
                        bubblecolor='blue'
                    else:
                        bubblecolor=bubblecolors[icase]
                    if bias[ivar,icase] > 0.0:
                        #linewidth=2.0
                        linewidth=0.0
                    else:
                        linewidth=0.0
                    circle = patches.Circle(
                            (x[icase],y[icase]),size[icase],
                            color=bubblecolor,
                            linewidth=linewidth,
                            alpha=0.20,
                        )
                    self.ax.add_patch(circle)

                    # draw the actual points
                    pointcolor=bubblecolor
                    circle = patches.Circle(
                            (x[icase],y[icase]),0.015,
                            color=pointcolor,
                            #edgecolor="black",
                            #linewidth=0.5,
                        )
                    self.ax.add_patch(circle)

                    # add labels to points if given
                    if varlabels is not None:
                        self.ax.text(
                                x[icase],y[icase]+labeloffset,
                                varlabels[ivar],
                                horizontalalignment='center',
                                verticalalignment='bottom',
                                color=bubblecolor,
                            )
        else:

            # transform coordinates
            r = ratio
            theta = numpy.arccos(cc)
            size = numpy.abs(bias)/2.0
            x,y = self.polar_transform(r,theta)
            for icase in range(len(x)):

                # draw "bubbles" around points
                if bubblecolors is None: 
                    bubblecolor='blue'
                else:
                    bubblecolor=bubblecolors[icase]
                if bias[icase] > 0.0:
                    #linewidth=2.0
                    linewidth=0.0
                else:
                    linewidth=0.0
                circle = patches.Circle(
                        (x[icase],y[icase]),size[icase],
                        color=bubblecolor,
                        linewidth=linewidth,
                        alpha=0.20,
                    )
                self.ax.add_patch(circle)

                # draw the actual points
                pointcolor=bubblecolor
                circle = patches.Circle(
                        (x[icase],y[icase]),0.015,
                        color=pointcolor,
                        #edgecolor="black",
                        #linewidth=0.5,
                    )
                self.ax.add_patch(circle)

                # add label to points if given
                if varlabels is not None:
                    self.ax.text(
                            x[icase],y[icase]+labeloffset,labels[icase],
                            horizontalalignment='center',
                            verticalalignment='bottom',
                            color=bubblecolor,
                        )

    def draw_rms(self,r,color='black',linestyle='dotted'):
        tmp = numpy.pi - numpy.arccos((r**2.0 + 1 - self.xymax**2.0)/2*r)
        theta1 = numpy.max([0.0,180.0/numpy.pi*tmp])
        theta2 = 180.0
        arc = patches.Arc(
                [1,0],2*r,2*r,
                theta1=theta1,
                theta2=theta2,
                color=color,
                linestyle=linestyle,
                linewidth=0.5,
            )
        self.ax.add_patch(arc)
