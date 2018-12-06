#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import h5py
from optparse import OptionParser
import glob
import requests
from PIL import Image
from io import BytesIO
#below is for an autorun test and can be deleted later...
import sys, os

def get_options():
    usage="%prog: [options]"

    parser = OptionParser(usage=usage)
    parser.add_option("-d", "--date", type="string", default="20140205",
                      help="Select date (YYYYMMDD) [default=%default]")
    parser.add_option("-t", "--time", type="string", default="1200",
                      help="Select time (hhmm) [default = %default]"
                      )
    parser.add_option("-f", "--prf", type=int, default=200,
                      help="Select pulse-repetition frequency (Hz) [default = %default]"
                      )
    (options, args) = parser.parse_args()

    return options

options = get_options()

time = float(options.time[0:2])*60 + float(options.time[2:4])
print time
time = int(time/15)
print time
time = time * 15
print time

hours = ['0']
minutes = ['0']
if (time / 60) < 10:
    hours.append(str(time / 60))
    hours = ''.join(hours)
else:
    hours = str(time / 60)
if (time % 60) < 10:
    minutes.append(str(time % 60))
    minutes = ''.join(minutes)
else:
    minutes = str(minutes % 60)
time = []
time.append(hours)
time.append(':')
time.append(minutes)
time.append(':00.000')
time = ''.join(time)
print time
print type(time)

payload1 = {'ursiCode':'EI764','year':options.date[0:4],'month':options.date[4:6],'day':options.date[6:8]}
print payload1
r1 = requests.get('http://car.uml.edu:8080/common/DIDBDayStationStatistic',params=payload1)
print r1.status_code
print r1.headers
print r1.encoding
print type(r1.text)
str1 = str(r1.text)
print type(str1)
time_idx = str1.find(time)
print time_idx
start_idx = str1.find(time)-120
print start_idx
search_str = str1[start_idx:time_idx]
print "search string is: ", search_str
start_idx = start_idx + search_str.find('mid=') + 4
print start_idx
end_idx = start_idx + 7
print end_idx
mid_str = str1[start_idx:end_idx]
print mid_str
payload2 = {'mid':mid_str}
r2 = requests.get('http://car.uml.edu:8080/common/ShowIonogram',params=payload2)
print r2.status_code
print r2.headers
print r2.encoding
print r2.url
print type(BytesIO(r2.content))
i = Image.open(BytesIO(r2.content))
print type(i)
i = np.asarray(i)
print type(i)

fstring = "ionogram." + options.date + "." + options.time + ".h5"
#print fstring
f = h5py.File("/home/pi/Desktop/1UAF_USRP-master/control_programs/swept_freq/"+fstring,'r')
#nfreqs = len(f.keys())/2;

dset = f['/OPower']
opower = np.empty([dset.shape[0], dset.shape[1]], float)
dset.read_direct(opower);

dset = f['/OVelocity']
ovelocity = np.empty([dset.shape[0], dset.shape[1]], float)
dset.read_direct(ovelocity);

dset = f['/XPower']
xpower = np.empty([dset.shape[0], dset.shape[1]], float)
dset.read_direct(xpower);

dset = f['/XVelocity']
xvelocity = np.empty([dset.shape[0], dset.shape[1]], float)
dset.read_direct(xvelocity);

freqs = f.attrs['Frequencies (kHz)']
minfreq = freqs[0]
maxfreq = freqs[freqs.shape[0]-1]
print "Start frequency:", minfreq
print "End frequency:", maxfreq
print "Number of frequencies:", freqs.shape[0]

ranges = f.attrs['Ranges (km)']
minrange = ranges[0]
print minrange
maxrange = ranges[ranges.shape[0]-1]
print "Start range:", minrange
print "End range:", maxrange
print "Number of range bins:", ranges.shape[0]

##ave_o = np.average(opower);
##ave_x = np.average(xpower);
##print ave_o
##print ave_x
##for i in range(0,nfreqs):
#	#opower[:,i] /= np.median(opower[:,i])
#	#xpower[:,i] /= np.median(xpower[:,i])
#	#image[:,i] /= np.average(image[:,i])
#
ext = [minfreq/1e3, maxfreq/1e3, minrange, maxrange]
ei_ext = [140,540,45,550]
asp = 2. * (maxfreq/1e3 - minfreq/1e3) / (maxrange - minrange);
#asp = 1
#print asp
#
opower = np.rot90(opower)
xpower = np.rot90(xpower)
ovelocity = np.rot90(ovelocity)
xvelocity = np.rot90(xvelocity)
##image[np.where(image < 0)] = 0
##image[np.where(image > 20)] = 20
ax1 = plt.subplot(121)
ax1.set_xlabel('Frequency (MHz)')
ax1.set_ylabel('Power (dB) at Range (km)')
title1 = ['PIonosonde Ionogram @ \n']
title1.append(options.date)
title1.append(' ')
title1.append(options.time)
title1 = ''.join(title1)
plt.title(title1)
plt.imshow(opower,extent=ext,aspect=asp,interpolation = "none")
plt.colorbar()

ax2 = plt.subplot(122)
ax2.set_xlabel('Frequency (MHz)')
title2 = ['Digisonde Ionogram @ \n']
title2.append(options.date)
title2.append(' ')
title2.append(time)
title2 = ''.join(title2)
plt.title(title2)
print type(i)
i = i[45:550,140:540,:]
plt.imshow(i)
plt.axis('off')
#plt.imshow(i)

plt.show()

#os.system("lxtask")
