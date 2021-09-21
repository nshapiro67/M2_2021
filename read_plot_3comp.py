#------------------------ importing basic packages
import matplotlib.pyplot as plt
import numpy as np

#------------------------ importing ObsPy functions
import obspy
from obspy import read
from obspy import UTCDateTime
from obspy.clients.fdsn import Client


plt.close("all")

evla = 16.982 
evlo = -99.773


#------------------------ selecting an FDSN datacenter
client = Client('IRIS')

#---------------------------------------------------- defining start time for the data
tstart = UTCDateTime("2021-09-08T01:47:47")

#-------------------- defining duration of the downloaded time series in sec
t_duration = 3*60*60

#--------------- network station
ntw = "IU"
sta = "CCM"
# ntw = "II"
# sta = "FFC"
# ntw = "IU"
# sta = "COLA"
# ntw = "G"
# sta = "SSB"

#--------------- downloading metadata
inv = client.get_stations(network=ntw,  starttime=tstart, endtime=tstart + t_duration, station=sta, channel="BH*", level="response")
#--------------- downloading waveforms
stin = client.get_waveforms(ntw, sta, "00", "BH*", tstart, tstart + t_duration, attach_response=True)

#--------- time exes series
dt = stin[0].stats.delta
npts = stin[0].stats.npts
time0 = dt*(np.linspace(1,npts,npts)-1)

dt = stin[1].stats.delta
npts = stin[1].stats.npts
time1 = dt*(np.linspace(1,npts,npts)-1)

dt = stin[2].stats.delta
npts = stin[2].stats.npts
time2 = dt*(np.linspace(1,npts,npts)-1)



#-----------------------  displacement seismograms 
pre_filt = (0.003, 0.005, 24.0, 25.0)                   # defining spectral band

stin.detrend('constant')
stin.detrend('linear')
stin.remove_response(output='DISP', pre_filt=pre_filt, inventory=inv)
stin.taper(.01)



#--------------- rotating seismograms to Z, N, and E components
stin.rotate(method="->ZNE", inventory=inv)


#------------------------------------------- plotting 3 displacement seismograms
fig, axs = plt.subplots(3, 1, sharex=True)
fig.subplots_adjust(hspace=0)

fig.canvas.set_window_title(stin[0].stats.starttime.ctime() + '     ' + stin[0].stats.station + '   displacement')

axs[0].plot(time0, stin[0].data)
axs[0].set_title(stin[0].stats.channel,loc='right', y=0.8)

axs[1].plot(time1, stin[1].data)
axs[1].set_title(stin[1].stats.channel, loc='right', y=0.8)
axs[1].set_ylabel("displacement (m)")

axs[2].plot(time2, stin[2].data)
axs[2].set_title(stin[2].stats.channel, loc='right', y=0.8)
axs[2].set_xlabel('time (s)')

plt.show()
#-------------------



#---------- rotation to NE
stla = inv[0][0].latitude
stlo = inv[0][0].longitude

(dst,baz,az) = obspy.geodetics.base.gps2dist_azimuth(stla,stlo,evla,evlo)

print('back azimuth ', baz, 'degrees')
print('distance ',dst/1000.,' km, (',dst/1000./111.11, 'degrees)')


stin.rotate(method="NE->RT", back_azimuth=baz)




#------------------------------------------- plotting 3 rotated displacement seismograms
fig, axs = plt.subplots(3, 1, sharex=True)
fig.subplots_adjust(hspace=0)

fig.canvas.set_window_title(stin[0].stats.starttime.ctime() + '     ' + stin[0].stats.station + '  rotated displacement')

axs[0].plot(time0, stin[0].data)
axs[0].set_title(stin[0].stats.channel,loc='right', y=0.8)

axs[1].plot(time1, stin[1].data)
axs[1].set_title(stin[1].stats.channel, loc='right', y=0.8)
axs[1].set_ylabel("displacement (m)")

axs[2].plot(time2, stin[2].data)
axs[2].set_title(stin[2].stats.channel, loc='right', y=0.8)
axs[2].set_xlabel('time (s)')

plt.show()
#-------------------


