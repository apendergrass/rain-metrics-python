#!/usr/bin/python
##########################################################################
### Angeline Pendergrass, January 18 2017.
### Starting from precipitation data, 
### 1. Calculate the distribution of rain
### 2. Plot the change from one climate state to another
### This code is ported from the matlab code shift-plus-increase-modes-demo, originally in matlab. 
### 
### You can read about these methods and cite the following papers about them: 
### Pendergrass, A.G. and D.L. Hartmann, 2014: Two modes of change of the                                                                                                       
###   distribution of rain. Journal of Climate, 27, 8357-8371.                                                                                                                  
###   doi:10.1175/JCLI-D-14-00182.1.                                                                                                                                            
### and the shift and increase modes of response of the rainfall distribution                                                                                                   
### to warming, occuring across ENSO events or global warming simulations.                                                                                                      
### The response to warming is described in:                                                                                                                                    
### Pendergrass, A.G. and D.L. Hartmann, 2014: Changes in the distribution                                                                                                      
###   of rain frequency and intensity in response to global warming.                                                                                                            
###   Journal of Climate, 27, 8372-8383. doi:10.1175/JCLI-D-14-00183.1.                                                                                                         
###
### See github.com/apendergrass for the latest info and updates. 
##########################################################################
import os
import sys
from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt

L=2.5e6 # % w/m2. latent heat of vaporization of water
wm2tommd=1./L*3600*24 # % conversion from w/m2 to mm/d

##### Load up the demo data. 
### It is daily average precipitation, in units of mm/d, with dimensions of lats, lons, and time. 
### From two time periods, as well as the change in global mean surface air temperature (a scalar) 
file1='pdistdemodata.nc' 
fh = Dataset(file1, mode='r')
lat = fh.variables['lat'][:]
lon = fh.variables['lon'][:]
pdata1 = fh.variables['pdata1'][:]
pdata2 = fh.variables['pdata2'][:]
dt=fh.variables['dt'][:]
fh.close()


def makedists(pdata,binl):
    ##### This is called from within makeraindist.
    ##### Caclulate distributions 
    pds=pdata.shape;    nlat=pds[1];    nlon=pds[0];    nd=pds[2]
    bins=np.append(0,binl)
    n=np.empty((nlon,nlat,len(binl)))
    binno=np.empty(pdata.shape)
    for ilon in range(nlon):
        for ilat in range(nlat):
            # this is the histogram - we'll get frequency from this
            thisn,thisbin=np.histogram(pdata[ilon,ilat,:],bins) 
            n[ilon,ilat,:]=thisn
            # these are the bin locations. we'll use these for the amount dist
            binno[ilon,ilat,:]=np.digitize(pdata[ilon,ilat,:],bins) - 1
    #### Calculate the number of days with non-missing data, for normalization
    ndmat=np.tile(np.expand_dims(np.nansum(n,axis=2),axis=2),(1,1,len(bins)-1))
    thisppdfmap=n/ndmat
    #### Iterate back over the bins and add up all the precip - this will be the rain amount distribution
    testpamtmap=np.empty(thisppdfmap.shape)
    for ibin in range(len(bins)-1):
        testpamtmap[:,:,ibin]=(pdata*(ibin==binno)).sum(axis=2)
    thispamtmap=testpamtmap/ndmat
    return thisppdfmap,thispamtmap


def makeraindist(pdata1,pdata2,lat,lon): 
    #### 1. Calculate bin structure. Note, these were chosen based on daily CMIP5 data - if you're doing something else you might want to change it
    sp1=pdata1.shape
    sp2=pdata2.shape
    if (sp1[1]!=len(lat))&(sp1[0]!=len(lon)):
        print 'pdata1 should be [days,lat,lon]'
    if (sp2[1]!=len(lat))&(sp2[0]!=len(lon)): 
        print 'pdata2 should be [days,lat,lon]'
    pmax=np.array([pdata1.max(),pdata2.max()]).max()/wm2tommd
    maxp=1500;# % choose an arbitrary upper bound for initial distribution, in w/m2
    minp=1;# % arbitrary lower bound, in w/m2. Make sure to set this low enough that you catch most of the rain. 
    #%%% thoughts: it might be better to specify the minimum threshold and the                                     
    #%%% bin spacing, which I have around 7%. The goals are to capture as much                                     
    #%%% of the distribution as possible and to balance sampling against                                           
    #%%% resolution. Capturing the upper end is easy: just extend the bins to                                      
    #%%% include the heaviest precipitation event in the dataset. The lower end                                    
    #%%% is harder: it can go all the way to machine epsilon, and there is no                                       
    #%%% obvious reasonable threshold for "rain" over a large spatial scale. The                                   
    #%%% value I chose here captures 97% of rainfall in CMIP5.                                                     
    nbins=100;
    binrlog=np.linspace(np.log(minp),np.log(maxp),nbins);
    dbinlog=np.diff(binrlog);
    binllog=binrlog-dbinlog[0];
    binr=np.exp(binrlog)/L*3600*24;
    binl=np.exp(binllog)/L*3600*24;
    dbin=dbinlog[0];
    binrlogex=binrlog;
    binrend=np.exp(binrlogex[len(binrlogex)-1])
    #% extend the bins until the maximum precip anywhere in the dataset falls
    #% within the bins
    # switch maxp to pmax if you want it to depend on your data
    while maxp>binr[len(binr)-1]:
        binrlogex=np.append(binrlogex,binrlogex[len(binrlogex)-1]+dbin)
        binrend=np.exp(binrlogex[len(binrlogex)-1]);
        binrlog=binrlogex;
        binllog=binrlog-dbinlog[0];
        binl=np.exp(binllog)/L*3600*24; #%% this is what we'll use to make distributions
        binr=np.exp(binrlog)/L*3600*24;
    bincrates=np.append(0,(binl+binr)/2)# % we'll use this for plotting.
    #### 2. Calculate distributions 
    ppdfmap,pamtmap=makedists(pdata1,binl);
    ppdfmap2,pamtmap2=makedists(pdata2,binl);
    #### 3. Spatially average distributions
    weight=np.tile(np.cos(lat*np.pi/180),(len(lon),1));
    weight=weight/weight.sum()
    weightp=np.tile(np.expand_dims(weight,axis=2),(1,1,ppdfmap.shape[2]))
    ppdf1=np.nansum(np.nansum(ppdfmap*weightp,axis=0),axis=0)
    pamt1=np.nansum(np.nansum(pamtmap*weightp,axis=0),axis=0)
    ppdf2=np.nansum(np.nansum(ppdfmap2*weightp,axis=0),axis=0)
    pamt2=np.nansum(np.nansum(pamtmap2*weightp,axis=0),axis=0)
    return ppdf1,pamt1,ppdf2,pamt2,bincrates


### Call the function to make the rain distribution
ppdf1, pamt1, ppdf2, pamt2, bincrates = makeraindist(pdata1,pdata2,lat,lon)

### Plot the change in rain amount and rain frequency distributions 
## you could turn this into a function if you wanted
#makedistplots(ppdf1,pamt1,ppdf2,pamt2,dt,bincrates)
#def makedistplots(ppdf1,pamt1,ppdf2,pamt2,dt,bincrates):

#### This is how we'll normalize to get changes per degree warming. 
pamt21k=pamt1+(pamt2-pamt1)/dt;
ppdf21k=ppdf1+(ppdf2-ppdf1)/dt;
ddry=((ppdf21k[0]-ppdf1[0])*100) # Change in dry days
    # % rain rates in mm/d for x axis ticks and labeling 
otn=np.linspace(1,9,9)
xtickrates=np.append(0,otn*.1)
xtickrates=np.append(xtickrates,otn)
xtickrates=np.append(xtickrates,otn*10)
xtickrates=np.append(xtickrates,otn*100)
xticks=np.interp(xtickrates,bincrates,range(0,len(bincrates))); #% bin numbers associated with nice number rain rate
xticks,indices=np.unique(xticks,return_index=True)
xtickrates=xtickrates[indices]
    ### Bin width - needed to normalize the rain amount distribution
db=(bincrates[2]-bincrates[1])/bincrates[1];
    ### Now we plot
plt.figure(figsize=(4,6))
plt.clf()
ax=plt.subplot(211)
plt.plot(range(0,len(pamt1)),(pamt21k-pamt1)/db, 'k')
plt.plot((0,len(pamt1)),(0,0),'0.5')
plt.ylim((-.05,.15))
plt.xlim((4,130))
    #plt.setp(ax,xticks=xticks,xticklabels=['0','0.1','','','','','','','','','','1','','','','','','','','','10','','','','','','','','','100','','','','','','','','','1000'])
plt.setp(ax,xticks=xticks,xticklabels=[''])
    #plt.xlabel('Rain rate (mm/d)')
plt.title('Rain amount change (mm/d/K)')
ax=plt.subplot(212)
plt.plot(range(0,len(ppdf1)),(ppdf21k-ppdf1)*100/db, 'k')
plt.plot((0,len(ppdf1)),(0,0),'0.5')
plt.ylim((-.06/db,.1/db))
plt.xlim((4,130))
    ### Annotate with the dry day frequency
t=plt.text(4,.095/db, "{:.1f}".format(ddry)+'%')
plt.setp(t,va='top',ha='left')
plt.setp(ax,xticks=xticks,xticklabels=['0','0.1','','','','','','','','','','1','','','','','','','','','10','','','','','','','','','100','','','','','','','','','1000'])
plt.xlabel('Rain rate (mm/d)')
plt.title('Rain frequency change (%/K)')
#plt.show()
filename="raindistdemo.pdf"
plt.savefig(filename)
print "wrote "+filename
plt.close()
#    return
