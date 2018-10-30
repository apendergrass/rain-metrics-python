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
### We'll just use the first time period here, but there are two, and also the change in global mean surface temperature. 
file1='pdistdemodata.nc' 
fh = Dataset(file1, mode='r')
lat = fh.variables['lat'][:]
lon = fh.variables['lon'][:]
pdata1 = fh.variables['pdata1'][:]
#pdata2 = fh.variables['pdata2'][:]
#dt=fh.variables['dt'][:]
fh.close()



def calcrainmetrics(pdistin,bincrates):
    ### This calculation can be applied to rain amount or rain frequency distributions 
    ### Here we'll do it for a distribution averaged over a region, but you could also do it at each grid point
    pdist=np.copy(pdistin)
    tile=np.array(0.1) # this is the threshold, 10% of rain amount or rain frequency
    pdist[0]=0 # If this is frequency, get rid of the dry frequency. If it's amount, it should already be zero or close to it.
    pmax=pdist.max()
    if pmax>0:
        imax=np.nonzero(pdist==pmax)
        rmax=np.interp(imax,range(0,len(bincrates)),bincrates)
        rainpeak=rmax[0][0]
        ### we're going to find the width by summing downward from pmax to lines at different heights, and then interpolating to figure out the rain rates that intersect the line. 
        theps=np.linspace(0.1,.99,99)*pmax
        thefrac=np.empty(theps.shape) 
        for i in range(len(theps)):
            thisp=theps[i]
            overp=(pdist-thisp)*(pdist> thisp);
            thefrac[i]=sum(overp)/sum(pdist)
        ptilerain=np.interp(-tile,-thefrac,theps) 
        #ptilerain/db ### check this against rain amount plot
        #ptilerain*100/db ### check this against rain frequency plot
        diffraintile=(pdist-ptilerain);
        alli=np.nonzero(diffraintile>0)
        afterfirst=alli[0][0]
        noistart=np.nonzero(diffraintile[0:afterfirst]<0)
        beforefirst=noistart[0][len(noistart[0])-1]
        incinds=range(beforefirst,afterfirst+1)
        ### need error handling on these for when inter doesn't behave well and there are multiple crossings
        if np.all(np.diff(diffraintile[incinds]) > 0):
            r1=np.interp(0,diffraintile[incinds],incinds) # this is ideally what happens. note: r1 is a bin index, not a rain rate. 
        else:
            r1=np.average(incinds) # in case interp won't return something meaningful, we use this kluge. 
        beforelast=alli[0][len(alli[0])-1]
        noiend=np.nonzero(diffraintile[beforelast:(len(diffraintile)-1)]<0)+beforelast
        afterlast=noiend[0][0]
        decinds=range(beforelast,afterlast+1)
        if np.all(np.diff(-diffraintile[decinds]) > 0):
            r2=np.interp(0,-diffraintile[decinds],decinds)
        else:
            r2=np.average(decinds) 
        ### Bin width - needed to normalize the rain amount distribution                                                                                                                
        db=(bincrates[2]-bincrates[1])/bincrates[1];
        rainwidth=(r2-r1)*db+1
        return rainpeak,rainwidth,(imax[0][0],pmax),(r1,r2,ptilerain)
    else:
        return 0,0,(0,pmax),(0,0,0)





    #### 1. Calculate bin structure. Note, these were chosen based on daily CMIP5 data - if you're doing something else you might want to change it
sp1=pdata1.shape
if (sp1[1]!=len(lat))&(sp1[0]!=len(lon)):
    print 'pdata1 should be [days,lat,lon]'
pmax=pdata1.max()/wm2tommd
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
            binno[ilon,ilat,:]=np.digitize(pdata[ilon,ilat,:],bins) 
    #### Calculate the number of days with non-missing data, for normalization
    ndmat=np.tile(np.expand_dims(np.nansum(n,axis=2),axis=2),(1,1,len(bins)-1))
    thisppdfmap=n/ndmat
    #### Iterate back over the bins and add up all the precip - this will be the rain amount distribution. 
    #### This step is probably the limiting factor and might be able to be made more efficient - I had a clever trick in matlab, but it doesn't work in python
    testpamtmap=np.empty(thisppdfmap.shape)
    for ibin in range(len(bins)-1):
        testpamtmap[:,:,ibin]=(pdata*(ibin==binno)).sum(axis=2)
    thispamtmap=testpamtmap/ndmat
    return thisppdfmap,thispamtmap

ppdfmap,pamtmap=makedists(pdata1,binl);

    #### 3. Spatially average distributions
weight=np.tile(np.cos(lat*np.pi/180),(len(lon),1));
weight=weight/weight.sum()
weightp=np.tile(np.expand_dims(weight,axis=2),(1,1,ppdfmap.shape[2]))
ppdf1=np.nansum(np.nansum(ppdfmap*weightp,axis=0),axis=0)
pamt1=np.nansum(np.nansum(pamtmap*weightp,axis=0),axis=0)



#### Calculate the rain metrics for the averaged distribution
rainamtpeak,rainamtwidth,plotpeakamt,plotwidthamt=calcrainmetrics(pamt1,bincrates)
rainpdfpeak,rainpdfwidth,plotpeakfreq,plotwidthfreq=calcrainmetrics(ppdf1,bincrates)
### Print the metrics for the averaged distribution
print "rain amount of average distribution"
print '  peak: '+"{:.1f}".format(rainamtpeak)+' mm/d'
print '  width: '+"{:.1f}".format(rainamtwidth)+' r2/r1'
print "rain frequency of average distribution"
print '  peak: '+"{:.1f}".format(rainpdfpeak)+' mm/d'
print '  width: '+"{:.1f}".format(rainpdfwidth)+' r2/r1'


#### Plot the average distribution and its metrics 

# This could be a function of its own
#makeshiftincplots(ppdf1,pamt1,ppdf2,pamt2,dt,bincrates);
#def makedistplots(ppdf1,pamt1,bincrates):

dry=ppdf1[0]*100 # Change in dry days
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
plt.plot(range(0,len(pamt1)),pamt1/db, 'k')
plt.plot(plotpeakamt[0],plotpeakamt[1]/db,'r*',markersize=10,mec='r')
plt.plot((plotwidthamt[0],plotwidthamt[1]),(plotwidthamt[2]/db,plotwidthamt[2]/db),'b');
plt.xlim((4,130))
plt.setp(ax,xticks=xticks,xticklabels=[''])
plt.title('Rain amount (mm/d)')
ax=plt.subplot(212)
plt.plot(range(0,len(ppdf1)),ppdf1/db*100, 'k')
plt.plot(plotpeakfreq[0],plotpeakfreq[1]/db*100,'r*',markersize=10,mec='r')
plt.plot((plotwidthfreq[0],plotwidthfreq[1]),(plotwidthfreq[2]*100/db,plotwidthfreq[2]*100/db),'b');
plt.xlim((4,130))
### Annotate with the dry day frequency
ymin, ymax = ax.get_ylim()
t=plt.text(4,ymax*0.95, "{:.1f}".format(dry)+'%')
plt.setp(t,va='top',ha='left')
plt.setp(ax,xticks=xticks,xticklabels=['0','0.1','','','','','','','','','','1','','','','','','','','','10','','','','','','','','','100','','','','','','','','','1000'])
plt.xlabel('Rain rate (mm/d)')
plt.title('Rain frequency (%)')
#plt.show()
filename="rainmetricdemo_averagedist.pdf"
plt.savefig(filename)
print "wrote "+filename
plt.close()

### now we'll make maps of the metrics to examine their spatial pattern. 

### Calculate the metrics for the distribution at each grid point
# ppdfmap,pamtmap=makedists(pdata1,binl);
amtpeakmap=np.empty((len(lon),len(lat)))
pdfpeakmap=np.empty((len(lon),len(lat)))
amtwidthmap=np.empty((len(lon),len(lat)))
pdfwidthmap=np.empty((len(lon),len(lat)))
for i in range(len(lon)):
    for j in range(len(lat)):
        rainpeak,rainwidth,plotpeak,plotwidth=calcrainmetrics(pamtmap[i,j,:],bincrates)
        amtpeakmap[i,j]=rainpeak
        amtwidthmap[i,j]=rainwidth
        rainpeak,rainwidth,plotpeak,plotwidth=calcrainmetrics(ppdfmap[i,j,:],bincrates)
        pdfpeakmap[i,j]=rainpeak
        pdfwidthmap[i,j]=rainwidth


## Now we'll plot them. The domain is the tropical north pacific. Note the abruptly high values of rain amount peak - those are Hawaii. 

def extents(f):
  delta = f[1] - f[0]
  return [f[0] - delta/2, f[-1] + delta/2]


plt.clf()
fig=plt.figure(figsize=(7,10))

ax=plt.subplot(321)
cax=plt.imshow(amtpeakmap.transpose(), interpolation='none',extent=extents(lon) + extents(lat), origin='lower',cmap='Greens')
cbar = plt.colorbar(cax,orientation='horizontal')
plt.title('Rain amount peak (mm/d)')

ax=plt.subplot(322)
cax=plt.imshow(pdfpeakmap.transpose(), interpolation='none',extent=extents(lon) + extents(lat), origin='lower',cmap='Greens')
cbar = plt.colorbar(cax,orientation='horizontal')
plt.title('Rain frequency peak (mm/d)')

ax=plt.subplot(323)
cax=plt.imshow(amtwidthmap.transpose(), interpolation='none',extent=extents(lon) + extents(lat), origin='lower',cmap='Greens')
cbar = plt.colorbar(cax,orientation='horizontal')
plt.title('Rain amount width (r2/r1)')

ax=plt.subplot(324)
cax=plt.imshow(pdfwidthmap.transpose(), interpolation='none',extent=extents(lon) + extents(lat), origin='lower',cmap='Greens')
cbar = plt.colorbar(cax,orientation='horizontal')
plt.title('Rain frequency width (r2/r1)')

ax=plt.subplot(325)
cax=plt.imshow(np.mean(pdata1,axis=2).transpose(), interpolation='none',extent=extents(lon) + extents(lat), origin='lower',cmap='Greens')
cbar = plt.colorbar(cax,orientation='horizontal')
plt.title('Average precip (mm/d)')

fig.text(0.5, 0.04, 'Longitude', ha='center')
fig.text(0.04, 0.5, 'Latitude', va='center', rotation='vertical')

#plt.show()                    



filename="rainmetricdemo_metricmaps.pdf"
plt.savefig(filename)
print "wrote "+filename
plt.close()
