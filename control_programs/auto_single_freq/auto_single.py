#!/usr/bin/env python
import time
from subprocess import call
from subprocess import check_output

def get_pid(name):
    return check_output(["pidof","s",name])

def get_ppid(pid):
    return check_output(["ps","--ppid",pid])

starttime=time.time()
filetime=time.localtime()
runcount = 0
while True:
    pid = get_pid("tilda")
    ppid = get_ppid(pid[0:4])
    #call(["renice", "-20", "-p", ppid])
    call(["renice", "-20", "-p", "2161"])
    timenow = time.localtime()
    filename = "ionogram.{:04d}{:02d}{:02d}.{:02d}{:02d}.h5".format(filetime.tm_year,filetime.tm_mon,filetime.tm_mday,filetime.tm_hour,filetime.tm_min)
    #call(["./single_freq", "--start", "5250", "--stop", "5350", "--npulses", "2", "--write", "--filename", filename, "--runcount", str(runcount)])
    call(["./single_freq", "--start", "3555", "--stop", "5000", "--npulses", "320", "--resolution", "5", "--first-range","50","--last-range", "750","--write"])
    runcount += 1
    time.sleep(180.0 - ((time.time() - starttime) % 180))
