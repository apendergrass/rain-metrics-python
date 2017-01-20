# rain-metrics-python
Python code to calculate the distribution of rain and metrics described in Pendergrass and Deser 2017

Starting from precipitation data,                                                                                                                                   
1. Calculate the distribution of rain                                                                                                                                     
2. Plot the change from one climate state to another                                                                           

Demo data are included in pdistdemodata.nc 

The rain distribution calculations are based on the Matlab shift-plus-increase-modes-demo.                              

You can read about the methods for calculating the distribution of rain here:

	Pendergrass, A.G. and D.L. Hartmann, 2014: Two modes of change of the 
  	distribution of rain. Journal of Climate, 27, 8357-8371. 
  	doi:10.1175/JCLI-D-14-00182.1.  
  
The response to warming is described in: 

	Pendergrass, A.G. and D.L. Hartmann, 2014: Changes in the distribution 
  	of rain frequency and intensity in response to global warming. 
  	Journal of Climate, 27, 8372-8383. doi:10.1175/JCLI-D-14-00183.1. 

The rain amount peak, rain frequency peak, rain amount width, and rain frequency width metrics are described in a forthcoming paper (currently being revised) for the Journal of Climate by Pendergrass and Deser. 

Please cite one or all of these papers. 

A bit of sample data from a member of the CESM1 large ensemble is included to get you going.
(https://www2.cesm.ucar.edu/models/experiments/LENS)

You can also use your own gridded precipiation dataset, of course.
For example, you can use daily rainfall data from a CMIP5 simulation
(GFDL-ESM2G is shown below, for the 1pctCO2 scenario)
which you might be able to download, for example from PCMDI
(https://pcmdi.llnl.gov/projects/cmip5/)

Or you can use data from GPCP 1dd
https://climatedataguide.ucar.edu/climate-data/gpcp-daily-global-precipitation-climatology-project
or TRMM 3B42
https://climatedataguide.ucar.edu/climate-data/trmm-tropical-rainfall-measuring-mission.

Get in touch if you have questions, or if you're interested in collaborating. 

*20 January 2017, Angeline Pendergrass, NCAR, Boulder CO. apgrass@ucar.edu*
