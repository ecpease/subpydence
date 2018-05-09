import flopy
import os
import numpy as np

mf = flopy.modflow.Modflow.load(os.path.join('HAGM.2012.nam'),model_ws='HAGM',load_only='dis')



dis = mf.dis



print(dis)


nstp = dis.nstp.array
perlen = dis.perlen.array

print(nstp)
print(perlen)


print(dis.nper)
filename = os.path.join('HAGM','HAGM.2012.hds')
headobj = flopy.utils.binaryfile.HeadFile(filename)
kstpker = np.array(headobj.get_kstpkper())
print(kstpker)


times = np.array(perlen)/365.25
print(times.cumsum()-10000)


newtimes = headobj.get_times()
newtimes = (np.array(newtimes)/365.25)-10000


print(newtimes)

