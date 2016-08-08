import numpy as numpy
import Scientific.IO.NetCDF as netcdf

class geo_var():

    def copy(self):
        newvar = map_var()
        newvar.data = numpy.ma.masked_array(self.data)
        newvar.data.mask = self.data.mask
        newvar._FillValue = self._FillValue
        newvar.lat = self.lat
        newvar.lon = self.lon
        return newvar

    def mask_area(self,cntl):
        self.data.mask = numpy.ma.mask_or(self.data.mask,cntl.data.mask)

    def lat_weights(self):
        wgt = numpy.cos(self.lat*numpy.pi/180.0)
        wgt = 2.0*wgt/wgt.sum()
        return wgt

    def weighted_area_average(self):
        wgt = self.lat_weights()
        ndims = len(self.data.shape)
        return numpy.ma.average(numpy.ma.average(self.data,ndims-2,wgt),ndims-2)

    def weighted_std(self):

        # get weights
        wgt_array = self.weight_array() 

        # make sure weight array is masked
        wgt_array = numpy.ma.masked_array(wgt_array)
        wgt_array.mask = self.data.mask

        # calculate standard deviation
        mean = self.weighted_area_average()
        num = numpy.ma.sum(wgt_array*(self.data-mean)**2.0)
        dem = numpy.ma.sum(wgt_array)
        return (num/dem)**0.5

    def weight_array(self):
        wgt_array = numpy.ma.zeros(self.data.shape)
        wgt_array.mask = self.data.mask
        wgt = self.lat_weights()
        for j in range(len(self.lat)):
            if len(self.data.shape) == 2:
                wgt_array[j,:] = wgt[j]
            elif len(self.data.shape) == 3:
                wgt_array[:,j,:] = wgt[j]
            elif len(self.data.shape) == 4:
                wgt_array[:,:,j,:] = wgt[j]
            elif len(self.data.shape) == 5:
                wgt_array[:,:,:,j,:] = wgt[j]
        return wgt_array

    def region_subset(self,(latmin,latmax),(lonmin,lonmax)):
        latarray = numpy.zeros(self.data.shape)
        lonarray = numpy.zeros(self.data.shape)
        if len(self.data.shape) == 2:
            for n in range(self.lat.size):
                latarray[n,:] = self.lat[n]
            for n in range(self.lon.size):
                lonarray[:,n] = self.lon[n]
        elif len(self.data.shape) == 3:
            for n in range(self.lat.size):
                latarray[:,n,:] = self.lat[n]
            for n in range(self.lon.size):
                lonarray[:,:,n] = self.lon[n]
        elif len(self.data.shape) == 4:
            for n in range(self.lat.size):
                latarray[:,:,n,:] = self.lat[n]
            for n in range(self.lon.size):
                lonarray[:,:,:,n] = self.lon[n]
        elif len(self.data.shape) == 5:
            for n in range(self.lat.size):
                latarray[:,:,:,n,:] = self.lat[n]
            for n in range(self.lon.size):
                lonarray[:,:,:,:,n] = self.lon[n]
        else:
            print "data shape not supported"
            exit()

        latmask = ~((latarray>=latmin) & (latarray<=latmax))
        lonmask = ~((lonarray>=lonmin) & (lonarray<=lonmax))

        self.data.mask = numpy.ma.mask_or(self.data.mask,latmask)
        self.data.mask = numpy.ma.mask_or(self.data.mask,lonmask)

    def regrid(self,latt,lont):
        """
        Method to regrid to a new target grid by only selecting those points
        in the source grid that have coordinates closest to those in the
        target grid. No smoothing is done by this method, so it should preserve
        the spatial correlation of the source data.
        """

        shape = self.data.shape
        if len(self.data.shape) == 2:
            data = numpy.ma.zeros([len(latt),len(lont)])
        elif len(self.data.shape) == 3:
            data = numpy.ma.zeros([shape[0],len(latt),len(lont)])
        elif len(self.data.shape) == 4:
            data = numpy.ma.zeros([shape[0],shape[1],len(latt),len(lont)])
        elif len(self.data.shape) == 5:
            data = numpy.ma.zeros([
                    shape[0],shape[1],shape[2],len(latt),len(lont)
                ])
        for i,lon in enumerate(lont):
            for j,lat in enumerate(latt):
                idxold = numpy.abs(self.lon-lon).argmin()
                idyold = numpy.abs(self.lat-lat).argmin()
                if len(shape) == 2:
                    data[j,i] = self.data[idyold,idxold]
                if len(shape) == 3:
                    data[:,j,i] = self.data[:,idyold,idxold]
                if len(shape) == 4:
                    data[:,:,j,i] = self.data[:,:,idyold,idxold]
                if len(shape) == 5:
                    data[:,:,:,j,i] = self.data[:,:,:,idyold,idxold]
        self.data = numpy.ma.masked_values(data,self._FillValue)
        self.lat = latt
        self.lon = lont


class map_var(geo_var):
    def __init__(self,filename=None,varname=None):
        if filename is not None and varname is not None:
            self.read(filename,varname)

    def read(self,filename,varname):

        # thick clouds
        if varname == "cltisccpthick":
            tmp = cltisccpthick_var(filename)
            self.data = numpy.ma.masked_values(tmp.data,tmp._FillValue)
            self._FillValue = tmp._FillValue
            self.lat = tmp.lat
            self.lon = tmp.lon
        elif varname == "cltmodisthick":
            tmp = cltmodisthick_var(filename)
            self.data = numpy.ma.masked_values(tmp.data,tmp._FillValue)
            self._FillValue = tmp._FillValue
            self.lat = tmp.lat
            self.lon = tmp.lon
        elif varname == "cltmisrthick":
            tmp = cltmisrthick_var(filename)
            self.data = numpy.ma.masked_values(tmp.data,tmp._FillValue)
            self._FillValue = tmp._FillValue
            self.lat = tmp.lat
            self.lon = tmp.lon
        elif varname == "cllisccpthick":
            tmp = cllisccpthick_var(filename)
            self.data = numpy.ma.masked_values(tmp.data,tmp._FillValue)
            self._FillValue = tmp._FillValue
            self.lat = tmp.lat
            self.lon = tmp.lon
        elif varname == "cllmodisthick":
            tmp = cllmodisthick_var(filename)
            self.data = numpy.ma.masked_values(tmp.data,tmp._FillValue)
            self._FillValue = tmp._FillValue
            self.lat = tmp.lat
            self.lon = tmp.lon
        elif varname == "cllmisrthick":
            tmp = cllmisrthick_var(filename)
            self.data = numpy.ma.masked_values(tmp.data,tmp._FillValue)
            self._FillValue = tmp._FillValue
            self.lat = tmp.lat
            self.lon = tmp.lon
        elif varname == "clmisccpthick":
            tmp = clmisccpthick_var(filename)
            self.data = numpy.ma.masked_values(tmp.data,tmp._FillValue)
            self._FillValue = tmp._FillValue
            self.lat = tmp.lat
            self.lon = tmp.lon
        elif varname == "clmmodisthick":
            tmp = clmmodisthick_var(filename)
            self.data = numpy.ma.masked_values(tmp.data,tmp._FillValue)
            self._FillValue = tmp._FillValue
            self.lat = tmp.lat
            self.lon = tmp.lon
        elif varname == "clmmisrthick":
            tmp = clmmisrthick_var(filename)
            self.data = numpy.ma.masked_values(tmp.data,tmp._FillValue)
            self._FillValue = tmp._FillValue
            self.lat = tmp.lat
            self.lon = tmp.lon
        elif varname == "clhisccpthick":
            tmp = clhisccpthick_var(filename)
            self.data = numpy.ma.masked_values(tmp.data,tmp._FillValue)
            self._FillValue = tmp._FillValue
            self.lat = tmp.lat
            self.lon = tmp.lon
        elif varname == "clhmodisthick":
            tmp = clhmodisthick_var(filename)
            self.data = numpy.ma.masked_values(tmp.data,tmp._FillValue)
            self._FillValue = tmp._FillValue
            self.lat = tmp.lat
            self.lon = tmp.lon
        elif varname == "clhmisrthick":
            tmp = clhmisrthick_var(filename)
            self.data = numpy.ma.masked_values(tmp.data,tmp._FillValue)
            self._FillValue = tmp._FillValue
            self.lat = tmp.lat
            self.lon = tmp.lon
        # other vars
        else:

            # open file
            ptr = netcdf.NetCDFFile(filename,"r")

            # read data
            data = ptr.variables[varname].getValue()[0,:,:]

            # read _FillValue, mask missing data
            self._FillValue = getattr(ptr.variables[varname],"_FillValue")
            self.data = numpy.ma.masked_values(data,self._FillValue)

            # read coordinate variables
            self.lat = ptr.variables["lat"].getValue()
            self.lon = ptr.variables["lon"].getValue()

class taylor_var(geo_var):
    def __init__(self,fileprefix,varname,cycle=True):
        if cycle:
            self.read_cycle(fileprefix,varname)
        else:
            self.read_data(fileprefix,varname)

    def read_data(self,filename,varname):

        # open file
        ptr = netcdf.NetCDFFile(filename,"r")

        # read coordinate variables
        self.lat = ptr.variables["lat"].getValue()
        self.lon = ptr.variables["lon"].getValue()

        # read data
        self.data = numpy.ma.zeros([1,self.lat.size,self.lon.size])
        tmp = map_var(filename,varname)
        self.data[0,:,:] = numpy.ma.masked_values(tmp.data,tmp._FillValue)
        self._FillValue = tmp._FillValue

    def read_cycle(self,fileprefix,varname):

        # open file to get dimension sizes
        filename = fileprefix+"01.nc"
        ptr = netcdf.NetCDFFile(filename,"r")
        self.lat = ptr.variables["lat"].getValue()
        self.lon = ptr.variables["lon"].getValue()

        # initialize array to hold data over the whole annual cycle
        self.data = numpy.ma.zeros([12,self.lat.size,self.lon.size])
        months = ["01","02","03","04","05","06","07","08","09","10","11","12"]

        # loop over months to populate data array
        for m,month in enumerate(months):

            # read data from single file
            filename = fileprefix+month+".nc"
            tmp = map_var(filename,varname)

            # populate array and mask missing values
            self._FillValue = tmp._FillValue
            self.data[m,:,:] = numpy.ma.masked_values(tmp.data,tmp._FillValue)

    def taylor_statistics(self,cntl):
        """Calculate Taylor statistics for making taylor diagrams."""

        # mask any missing values
        self.data.mask = numpy.ma.mask_or(self.data.mask,cntl.data.mask)
        cntl.data.mask = numpy.ma.mask_or(cntl.data.mask,self.data.mask)

        # copy weights to full data-dimensioned array
        wgt = self.weight_array()

        # make sure weight array is masked
        wgt.mask = numpy.ma.mask_or(self.data.mask,cntl.data.mask)

        # calculate sums and means
        sumwgt = numpy.ma.sum(wgt)
        meanself = numpy.ma.sum(wgt*self.data)/sumwgt
        meancntl = numpy.ma.sum(wgt*cntl.data)/sumwgt

        # calculate variances
        stdself = (numpy.ma.sum(wgt*(self.data-meanself)**2.0)/sumwgt)**0.5
        stdcntl = (numpy.ma.sum(wgt*(cntl.data-meancntl)**2.0)/sumwgt)**0.5

        # calculate correlation coefficient
        ccnum = numpy.ma.sum(wgt*(self.data-meanself)*(cntl.data-meancntl))
        ccdem = sumwgt*stdself*stdcntl
        cc = ccnum/ccdem

        # calculate variance ratio
        ratio = stdself/stdcntl

        # calculate bias
        bias = (meanself - meancntl)/numpy.abs(meancntl)

        # calculate centered pattern RMS difference
        rmssum = numpy.ma.sum(wgt*((self.data-meanself)-(cntl.data-meancntl))**2.0)
        rmserr = (rmssum/sumwgt)**0.5
        rmsnorm = rmserr/stdcntl

        return ratio,cc,bias,rmsnorm

class hist_var(geo_var):
    def __init__(self,filename,varname,time=0):
        self.ptr = netcdf.NetCDFFile(filename,"r")
        self.read(varname,time=time)

    def read(self,varname,time=0):
        # read data
        data = self.ptr.variables[varname].getValue()[time,:,:,:,:]
        self._FillValue = getattr(self.ptr.variables[varname],"_FillValue")
        self.data = numpy.ma.masked_values(data,self._FillValue)

        # read coordinate variables
        self.lat = self.ptr.variables["lat"].getValue()
        self.lon = self.ptr.variables["lon"].getValue()
        if varname == "clMISR":
            self.tau = self.ptr.variables["tau"].getValue()
            self.cth = self.ptr.variables["cth"].getValue()
        elif varname == "clisccp":
            self.tau = self.ptr.variables["tau"].getValue()
            self.plev = self.ptr.variables["plev"].getValue()
        elif varname == "clmodis":
            self.tau = self.ptr.variables["tau"].getValue()
            self.plev = self.ptr.variables["plev"].getValue()

    def mask_area(self,cntl):
        # NOTE: this is ugly, and I don't think it is working right...
        # maybe try to interpolate mask data to self data grid and then
        # use a mask_or like in map_var
        area_mask = cntl.data.mask[0,0,:,:]
        for i in range(self.data.shape[0]):
            for j in range(self.data.shape[1]):
                self.data.mask[i,j,:,:] = (area_mask[:,:] | self.data.mask[i,j,:,:])


class tau_var(geo_var):
    def __init__(self,filename,varname):
        self.ptr = netcdf.NetCDFFile(filename,"r")
        self.read(varname)

    def read(self,varname):
        data = self.ptr.variables[varname].getValue()[0,:,:,:,:]
        self._FillValue = getattr(self.ptr.variables[varname],"_FillValue")
        self.data = numpy.ma.masked_values(data,self._FillValue)
        self.varname = varname

        # coordinate variables
        self.tau = self.ptr.variables["tau"].getValue()
        self.lat = self.ptr.variables["lat"].getValue()
        self.lon = self.ptr.variables["lon"].getValue()
        if "cth" in self.ptr.variables.keys():
            self.cth = self.ptr.variables["cth"].getValue()
        if "plev" in self.ptr.variables.keys():
            self.plev = self.ptr.variables["plev"].getValue()

    def mask_area(self,cntl):
        mask = numpy.array(self.data<0.0)
        for idx in range(mask.shape[0]):
            for idy in range(mask.shape[1]):
                mask[idx,idy,:,:] = cntl.data.mask[0,0,:,:]
        self.data.mask = numpy.ma.mask_or(self.data.mask,mask)

    def calculate(self):
        lat = self.ptr.variables["lat"].getValue()
        wgt = self.lat_weights()
        matrix = numpy.ma.average(numpy.ma.average(self.data,2,wgt),2)

        tau = self.ptr.variables["tau"].getValue()
        imin = next(idx for idx,value in enumerate(tau) if value > 1.3)
        if self.varname == "clMISR":
            self.histogram = matrix.sum(1)[imin:]
        else:
            self.histogram = matrix.sum(0)[imin:]

class cth_var(geo_var):
    def __init__(self,filename,varname):
        self.ptr = netcdf.NetCDFFile(filename,"r")
        self.read(varname)
        # self.calculate()

    def read(self,varname):
        # full joint histogram data
        data = self.ptr.variables[varname].getValue()[0,:,:,:,:]
        self._FillValue = getattr(self.ptr.variables[varname],"_FillValue")
        self.data = numpy.ma.masked_values(data,self._FillValue)

        # coordinate varibles
        self.cth = self.ptr.variables["cth"].getValue()
        self.tau = self.ptr.variables["tau"].getValue()
        self.lat = self.ptr.variables["lat"].getValue()
        self.lon = self.ptr.variables["lon"].getValue()
        self.cth_bnds = self.ptr.variables["cth_bnds"].getValue()[1:,:]


    def mask_area(self,cntl):
        mask = numpy.array(self.data<0.0)
        for idx in range(mask.shape[0]):
            for idy in range(mask.shape[1]):
                mask[idx,idy,:,:] = cntl.data.mask[0,0,:,:] 
        self.data.mask = numpy.ma.mask_or(self.data.mask,mask)

    def calculate(self):
        wgt = self.lat_weights()
        matrix = numpy.ma.average(numpy.ma.average(self.data,2,wgt),2)

        cth = self.ptr.variables["cth"].getValue()
        tau = self.ptr.variables["tau"].getValue()
        imin = next(idx for idx,value in enumerate(tau) if value > 1.3)
        self.histogram = matrix[imin:,:].sum(0)[1:]
        self.combine_bins()

    def combine_bins(self):
        # combine bins
        newcthbnds = numpy.zeros([self.cth_bnds.shape[0]-2,2])
        newhistogram = numpy.zeros([self.histogram.shape[0]-2])
        newhistogram[0] = self.histogram[0]
        newhistogram[1] = self.histogram[1]
        newhistogram[2] = self.histogram[2]+self.histogram[3]
        newhistogram[3] = self.histogram[4]+self.histogram[5]
        newhistogram[4:] = self.histogram[6:]
        self.histogram = newhistogram

        newcthbnds[0,:] = self.cth_bnds[0,:]
        newcthbnds[1,:] = self.cth_bnds[1,:]
        newcthbnds[2,0] = self.cth_bnds[2,0]
        newcthbnds[2,1] = self.cth_bnds[3,1]
        newcthbnds[3,0] = newcthbnds[2,1]
        newcthbnds[3,1] = self.cth_bnds[5,1]
        newcthbnds[4:,:] = self.cth_bnds[6:,:]
        self.cth_bnds = newcthbnds

class clcth_var(geo_var):
    def __init__(self,filename,varname):
        self.ptr = netcdf.NetCDFFile(filename,"r")
        self.read(varname)
        # self.calculate()

    def read(self,varname):
        # full joint histogram data
        data = self.ptr.variables[varname].getValue()[0,:,:,:,:]
        self._FillValue = getattr(self.ptr.variables[varname],"_FillValue")
        self.data = numpy.ma.masked_values(data,self._FillValue)

        # cloud top height bounds
        self.cth_bnds = self.ptr.variables["cth_bnds"].getValue()
        self.lat = self.ptr.variables["lat"].getValue()
        self.lon = self.ptr.variables["lon"].getValue()
        self.tau = self.ptr.variables["tau"].getValue()
        self.cth = self.ptr.variables["cth"].getValue()

    def mask_area(self,cntl):
        mask = numpy.array(self.data<0.0)
        for idx in range(mask.shape[0]):
            for idy in range(mask.shape[1]):
                mask[idx,idy,:,:] = cntl.data.mask[0,0,:,:] 
        self.data.mask = numpy.ma.mask_or(self.data.mask,mask)

    def calculate(self):
        tminind = next(idx for idx,value in enumerate(self.tau) if value > 1.3)
        hminind = next(idx for idx,value in enumerate(self.cth) if value > 0.0)
        self.histogram = numpy.ma.sum(self.data[tminind:,hminind:,:,:],0)
        self.tau = self.tau[tminind:]
        self.cth = self.cth[hminind:]
        self.cth_bnds = self.cth_bnds[hminind:,:]

    def combine_bins(self):
        # combine bins
        newcthbnds = numpy.zeros([self.cth_bnds.shape[0]-2,2])
        newhistogram = numpy.zeros([
                self.histogram.shape[0]-2,
                self.histogram.shape[1],
                self.histogram.shape[2]
            ])
        newhistogram[0,:,:] = self.histogram[0,:,:]
        newhistogram[1,:,:] = self.histogram[1,:,:]
        newhistogram[2,:,:] = self.histogram[2,:,:]+self.histogram[3,:,:]
        newhistogram[3,:,:] = self.histogram[4,:,:]+self.histogram[5,:,:]
        newhistogram[4:,:,:] = self.histogram[6:,:,:]
        self.histogram = newhistogram

        newcthbnds[0,:] = self.cth_bnds[0,:]
        newcthbnds[1,:] = self.cth_bnds[1,:]
        newcthbnds[2,0] = self.cth_bnds[2,0]
        newcthbnds[2,1] = self.cth_bnds[3,1]
        newcthbnds[3,0] = newcthbnds[2,1]
        newcthbnds[3,1] = self.cth_bnds[5,1]
        newcthbnds[4:,:] = self.cth_bnds[6:,:]
        self.cth_bnds = newcthbnds

class cl_var(geo_var):
    def __init__(self,filename,varname):
        self.varname = varname
        self.ptr = netcdf.NetCDFFile(filename,"r")
        self.read(varname)

    def read(self,varname):
        data = self.ptr.variables[varname].getValue()[0,:,:,:]
        if "_FillValue" in self.ptr.variables.keys():
            self._FillValue = getattr(self.ptr.variables[varname],"_FillValue")
            self.data = numpy.ma.masked_values(data,self._FillValue)
        else:
            self.data = data

        # coordinate variables
        self.lev = self.ptr.variables["lev"].getValue()
        self.lat = self.ptr.variables["lat"].getValue()
        self.lon = self.ptr.variables["lon"].getValue()

class cltisccp_var(geo_var):
    def __init__(self,filename):
        clisccp = hist_var(filename,"clisccp")
        tind1 = next(idx for idx,val in enumerate(clisccp.tau) if val > 0.3)
        self.data = numpy.ma.sum(numpy.ma.sum(clisccp.data[:,tind1:,:,:],1),0)
        self.lat = clisccp.lat
        self.lon = clisccp.lon
        self._FillValue = clisccp._FillValue


class cltmodis_var(geo_var):
    def __init__(self,filename):
        clmodis = hist_var(filename,"clmodis")
        tind1 = next(idx for idx,val in enumerate(clmodis.tau) if val > 0.3)
        self.data = numpy.ma.sum(numpy.ma.sum(clmodis.data[:,tind1:,:,:],1),0)
        self.lat = clmodis.lat
        self.lon = clmodis.lon
        self._FillValue = clmodis._FillValue


class cltmisr_var(geo_var):
    def __init__(self,filename):
        clMISR = hist_var(filename,"clMISR")
        tind1 = next(idx for idx,val in enumerate(clMISR.tau) if val > 0.3)
        self.data = numpy.ma.sum(numpy.ma.sum(clMISR.data[tind1:,:,:,:],1),0)
        self.lat = clMISR.lat
        self.lon = clMISR.lon
        self._FillValue = clMISR._FillValue


class cllisccp_var(geo_var):
    def __init__(self,filename):
        clisccp = hist_var(filename,"clisccp")
        t1 = next(idx for idx,val in enumerate(clisccp.tau) if val > 0.3)
        p2 = next(idx for idx,val in enumerate(clisccp.plev) if val < 68000)
        self.data = numpy.ma.sum(numpy.ma.sum(clisccp.data[:p2,t1:,:,:],1),0)
        self.lat = clisccp.lat
        self.lon = clisccp.lon
        self._FillValue = clisccp._FillValue


class cllmodis_var(geo_var):
    def __init__(self,filename):
        clmodis = hist_var(filename,"clmodis")
        t1 = next(idx for idx,val in enumerate(clmodis.tau) if val > 0.3)
        p2 = next(idx for idx,val in enumerate(clmodis.plev) if val < 68000)
        self.data = numpy.ma.sum(numpy.ma.sum(clmodis.data[:p2,t1:,:,:],1),0)
        self.lat = clmodis.lat
        self.lon = clmodis.lon
        self._FillValue = clmodis._FillValue


class cllmisr_var(geo_var):
    def __init__(self,filename):
        clMISR = hist_var(filename,"clMISR")
        t1 = next(idx for idx,val in enumerate(clMISR.tau) if val > 0.3)
        p1 = next(idx for idx,val in enumerate(clMISR.cth) if val > 0.0)
        p2 = next(idx for idx,val in enumerate(clMISR.cth) if val > 3000.0)
        self.data = numpy.ma.sum(numpy.ma.sum(clMISR.data[t1:,p1:p2,:,:],1),0)
        self.lat = clMISR.lat
        self.lon = clMISR.lon
        self._FillValue = clMISR._FillValue


class clmisccp_var(geo_var):
    def __init__(self,filename):
        clisccp = hist_var(filename,"clisccp")
        t1 = next(idx for idx,val in enumerate(clisccp.tau) if val > 0.3)
        p1 = next(idx for idx,val in enumerate(clisccp.plev) if val < 68000)
        p2 = next(idx for idx,val in enumerate(clisccp.plev) if val < 44000)
        self.data = numpy.ma.sum(
                numpy.ma.sum(clisccp.data[p1:p2,t1:,:,:],1),0
            )
        self.lat = clisccp.lat
        self.lon = clisccp.lon
        self._FillValue = clisccp._FillValue


class clmmodis_var(geo_var):
    def __init__(self,filename):
        clmodis = hist_var(filename,"clmodis")
        t1 = next(idx for idx,val in enumerate(clmodis.tau) if val > 0.3)
        p1 = next(idx for idx,val in enumerate(clmodis.plev) if val < 68000)
        p2 = next(idx for idx,val in enumerate(clmodis.plev) if val < 44000)
        self.data = numpy.ma.sum(
                numpy.ma.sum(clmodis.data[p1:p2,t1:,:,:],1),0
            )
        self.lat = clmodis.lat
        self.lon = clmodis.lon
        self._FillValue = clmodis._FillValue


class clmmisr_var(geo_var):
    def __init__(self,filename):
        clMISR = hist_var(filename,"clMISR")
        t1 = next(idx for idx,val in enumerate(clMISR.tau) if val > 0.3)
        p1 = next(idx for idx,val in enumerate(clMISR.cth) if val > 3000.0)
        p2 = next(idx for idx,val in enumerate(clMISR.cth) if val > 7000.0)
        self.data = numpy.ma.sum(numpy.ma.sum(clMISR.data[t1:,p1:p2,:,:],1),0)
        self.lat = clMISR.lat
        self.lon = clMISR.lon
        self._FillValue = clMISR._FillValue


class clhisccp_var(geo_var):
    def __init__(self,filename):
        clisccp = hist_var(filename,"clisccp")
        t1 = next(idx for idx,val in enumerate(clisccp.tau) if val > 0.3)
        p1 = next(idx for idx,val in enumerate(clisccp.plev) if val < 44000)
        self.data = numpy.ma.sum(
                numpy.ma.sum(clisccp.data[p1:,t1:,:,:],1),0
            )
        self.lat = clisccp.lat
        self.lon = clisccp.lon
        self._FillValue = clisccp._FillValue


class clhmodis_var(geo_var):
    def __init__(self,filename):
        clmodis = hist_var(filename,"clmodis")
        t1 = next(idx for idx,val in enumerate(clmodis.tau) if val > 0.3)
        p1 = next(idx for idx,val in enumerate(clmodis.plev) if val < 44000)
        self.data = numpy.ma.sum(
                numpy.ma.sum(clmodis.data[p1:,t1:,:,:],1),0
            )
        self.lat = clmodis.lat
        self.lon = clmodis.lon
        self._FillValue = clmodis._FillValue


class clhmisr_var(geo_var):
    def __init__(self,filename):
        clMISR = hist_var(filename,"clMISR")
        t1 = next(idx for idx,val in enumerate(clMISR.tau) if val > 0.3)
        p1 = next(idx for idx,val in enumerate(clMISR.cth) if val > 7000.0)
        self.data = numpy.ma.sum(numpy.ma.sum(clMISR.data[t1:,p1:,:,:],1),0)
        self.lat = clMISR.lat
        self.lon = clMISR.lon
        self._FillValue = clMISR._FillValue


class cltisccpthick_var(geo_var):
    def __init__(self,filename):
        clisccp = hist_var(filename,"clisccp")
        t1 = next(idx for idx,val in enumerate(clisccp.tau) if val > 23.0)
        self.data = numpy.ma.sum(numpy.ma.sum(clisccp.data[:,t1:,:,:],1),0)
        self.lat = clisccp.lat
        self.lon = clisccp.lon
        self._FillValue = clisccp._FillValue


class cltmodisthick_var(geo_var):
    def __init__(self,filename):
        clmodis = hist_var(filename,"clmodis")
        t1 = next(idx for idx,val in enumerate(clmodis.tau) if val > 23.0)
        self.data = numpy.ma.sum(numpy.ma.sum(clmodis.data[:,t1:,:,:],1),0)
        self.lat = clmodis.lat
        self.lon = clmodis.lon
        self._FillValue = clisccp._FillValue


class cltmisrthick_var(geo_var):
    def __init__(self,filename):
        clMISR = hist_var(filename,"clMISR")
        t1 = next(idx for idx,val in enumerate(clMISR.tau) if val > 23.0)
        self.data = numpy.ma.sum(numpy.ma.sum(clMISR.data[t1:,:,:,:],1),0)
        self.lat = clMISR.lat
        self.lon = clMISR.lon
        self._FillValue = clMISR._FillValue


class cllisccpthick_var(geo_var):
    def __init__(self,filename):
        clisccp = hist_var(filename,"clisccp")
        t1 = next(idx for idx,val in enumerate(clisccp.tau) if val > 23.0)
        p2 = next(idx for idx,val in enumerate(clisccp.plev) if val < 68000)
        self.data = numpy.ma.sum(numpy.ma.sum(clisccp.data[:p2,t1:,:,:],1),0)
        self.lat = clisccp.lat
        self.lon = clisccp.lon
        self._FillValue = clisccp._FillValue


class cllmodisthick_var(geo_var):
    def __init__(self,filename):
        clmodis = hist_var(filename,"clmodis")
        t1 = next(idx for idx,val in enumerate(clmodis.tau) if val > 23.0)
        p2 = next(idx for idx,val in enumerate(clmodis.plev) if val < 68000)
        self.data = numpy.ma.sum(numpy.ma.sum(clmodis.data[:p2,t1:,:,:],1),0)
        self.lat = clmodis.lat
        self.lon = clmodis.lon
        self._FillValue = clmodis._FillValue


class cllmisrthick_var(geo_var):
    def __init__(self,filename):
        clMISR = hist_var(filename,"clMISR")
        t1 = next(idx for idx,val in enumerate(clMISR.tau) if val > 23.0)
        p1 = next(idx for idx,val in enumerate(clMISR.cth) if val > 0.0)
        p2 = next(idx for idx,val in enumerate(clMISR.cth) if val > 3000.0)
        self.data = numpy.ma.sum(numpy.ma.sum(clMISR.data[t1:,p1:p2,:,:],1),0)
        self.lat = clMISR.lat
        self.lon = clMISR.lon
        self._FillValue = clMISR._FillValue


class clmisccpthick_var(geo_var):
    def __init__(self,filename):
        clisccp = hist_var(filename,"clisccp")
        t1 = next(idx for idx,val in enumerate(clisccp.tau) if val > 23.0)
        p1 = next(idx for idx,val in enumerate(clisccp.plev) if val < 68000)
        p2 = next(idx for idx,val in enumerate(clisccp.plev) if val < 44000)
        self.data = numpy.ma.sum(
                numpy.ma.sum(clisccp.data[p1:p2,t1:,:,:],1),0
            )
        self.lat = clisccp.lat
        self.lon = clisccp.lon
        self._FillValue = clisccp._FillValue


class clmmodisthick_var(geo_var):
    def __init__(self,filename):
        clmodis = hist_var(filename,"clmodis")
        t1 = next(idx for idx,val in enumerate(clmodis.tau) if val > 23.0)
        p1 = next(idx for idx,val in enumerate(clmodis.plev) if val < 68000)
        p2 = next(idx for idx,val in enumerate(clmodis.plev) if val < 44000)
        self.data = numpy.ma.sum(
                numpy.ma.sum(clmodis.data[p1:p2,t1:,:,:],1),0
            )
        self.lat = clmodis.lat
        self.lon = clmodis.lon
        self._FillValue = clmodis._FillValue


class clmmisrthick_var(geo_var):
    def __init__(self,filename):
        clMISR = hist_var(filename,"clMISR")
        t1 = next(idx for idx,val in enumerate(clMISR.tau) if val > 23.0)
        p1 = next(idx for idx,val in enumerate(clMISR.cth) if val > 3000.0)
        p2 = next(idx for idx,val in enumerate(clMISR.cth) if val > 7000.0)
        self.data = numpy.ma.sum(numpy.ma.sum(clMISR.data[t1:,p1:p2,:,:],1),0)
        self.lat = clMISR.lat
        self.lon = clMISR.lon
        self._FillValue = clMISR._FillValue


class clhisccpthick_var(geo_var):
    def __init__(self,filename):
        clisccp = hist_var(filename,"clisccp")
        t1 = next(idx for idx,val in enumerate(clisccp.tau) if val > 23.0)
        p1 = next(idx for idx,val in enumerate(clisccp.plev) if val < 44000)
        self.data = numpy.ma.sum(
                numpy.ma.sum(clisccp.data[p1:,t1:,:,:],1),0
            )
        self.lat = clisccp.lat
        self.lon = clisccp.lon
        self._FillValue = clisccp._FillValue


class clhmodisthick_var(geo_var):
    def __init__(self,filename):
        clmodis = hist_var(filename,"clmodis")
        t1 = next(idx for idx,val in enumerate(clmodis.tau) if val > 23.0)
        p1 = next(idx for idx,val in enumerate(clmodis.plev) if val < 44000)
        self.data = numpy.ma.sum(
                numpy.ma.sum(clmodis.data[p1:,t1:,:,:],1),0
            )
        self.lat = clmodis.lat
        self.lon = clmodis.lon
        self._FillValue = clmodis._FillValue


class clhmisrthick_var(geo_var):
    def __init__(self,filename):
        clMISR = hist_var(filename,"clMISR")
        t1 = next(idx for idx,val in enumerate(clMISR.tau) if val > 23.0)
        p1 = next(idx for idx,val in enumerate(clMISR.cth) if val > 7000.0)
        self.data = numpy.ma.sum(numpy.ma.sum(clMISR.data[t1:,p1:,:,:],1),0)
        self.lat = clMISR.lat
        self.lon = clMISR.lon
        self._FillValue = clMISR._FillValue
