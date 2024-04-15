import numpy as np

topofile = 'tmp2/teton.asc'
topofile_new = 'TetonLarge.topo'


f = open(topofile,'r')
fnew = open(topofile_new,'w')

l = f.readline()
ncols = np.fromstring(l.split()[1].strip(),sep=' ')                    
fnew.write("%s %s\n" % (l.split()[1],l.split()[0]))

l = f.readline()
nrows = np.fromstring(l.split()[1].strip(),sep=' ')
fnew.write("%s %s\n" % (l.split()[1],l.split()[0]))

l = f.readline()
xllcorner = np.fromstring(l.split()[1].strip(),sep=' ')
fnew.write("%s %s\n" % (l.split()[1],l.split()[0]))

l = f.readline()
yllcorner = np.fromstring(l.split()[1].strip(),sep=' ')
fnew.write("%s %s\n" % (l.split()[1],l.split()[0]))

l = f.readline()
cellsize = np.fromstring(l.split()[1].strip(),sep=' ')
fnew.write("%s %s\n" % (l.split()[1],l.split()[0]))

l = f.readline()
NODATA_value  = np.fromstring(l.split()[1].strip(),sep=' ')
fnew.write("%s %s\n" % (l.split()[1],l.split()[0]))

# This loop writes out data as one long column, which I think GeoClaw needs.
k = 0
for j in range(0,nrows):
    l = f.readline()
    datastr  = l.split()
    for i in range(0,ncols):
        fnew.write("%s\n" % datastr[i])
        k = k + 1

fnew.close()

print("ncols = %8d" % ncols)
print("nrows = %8d" % nrows)
print("xll = %20.12f" % xllcorner)
print("yll = %20.12f" % yllcorner)
print("cellsize = %20.12f" % cellsize)
print("nodata = %g" % NODATA_value)
print("") 
print("Data elements written out : %d" % k)
print("Data elements expected    : %d" % (nrows*ncols))


f.close()


