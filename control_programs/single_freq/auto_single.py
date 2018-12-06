import time
from subprocess import call
starttime=time.time()
while True:
    call(["./2single_freq", "--start", "5000", "--stop", "5000", "--write"])
    time.sleep(60.0 - ((time.time() - starttime) % 60))
