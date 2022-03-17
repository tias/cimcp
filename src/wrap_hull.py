#!/usr/bin/python -u
# run FIM_CP
#
# ./run_fimcp.py
# Tias, 2008-02-26
import os
import sys
import math
# If the package has been installed correctly, this should work:
from numpy import *
import Gnuplot, Gnuplot.funcutils

def getname(file):
    return os.path.splitext(os.path.split(file)[1])[0]

BINARYDIR='./'
DATADIR='../uci_d8b/'
name=getname(sys.argv[0])

if len(sys.argv) == 1:
    print """Usage: %s <datafile>
\tdatafile is a fimi-style dataset with the class label in the last column"""%(sys.argv[0])
    exit(1)

# read datafile name
datafile = sys.argv[1]
if (not os.path.isfile(datafile)):
    print "Error, \"%s\" is not a file"%(datafile)
    exit(1)


#### Running model ####
binary = './'+BINARYDIR+'cimcp_convexhull'
cmd = "%s -datafile %s -output normal"%(binary, datafile)
print 'Running:\t'+cmd
(_,fp) = os.popen4(cmd, 'r')
data = []
loading = False
for line in fp:
    if (line == "Convex Hull:\n"):
        data = []
        loading = True
        continue
    if loading:
        if line[0] != "\t":
            loading = False
            continue
        data.append(line.strip())

fp.close()
#print data


#### Calculating measures, plotting data ####
# calculate chi2
def chi2(X, a, Y, b):
    yaminxb = Y*a - X*b
    check = float((a+b)*(X+Y-a-b))
    if (check == 0):
        #print "warning, devision by 0 caught"
        return 0
    one = yaminxb / check
    two = yaminxb / float(X*Y)
    calc = one*two*(X+Y)
    
    return calc
# calculate infgain
def nlogn(n):
    if (n == 0):
        return 0
    return -math.log(n)*n
def infgain(X, pos_x, Y, pos_y):
    tot = float(X+Y)
    base = nlogn(X/tot) + nlogn(Y/tot)

    pos = float(pos_x + pos_y)
    neg = tot - pos
    calcPos = 0
    if (pos != 0):
        calcPos = pos * ( nlogn(    pos_x/pos) + nlogn(    pos_y/pos) )
    calcNeg = 0
    if (neg != 0):
        calcNeg = neg * ( nlogn((X-pos_x)/neg) + nlogn((Y-pos_y)/neg) )
    val = base - (calcPos+calcNeg)/tot
    #print "%.5f - (%.5f/%.5f) = %.5f"%(base,calc,tot,val)
    return val

ULpoints = []
LRpoints = []
inUL = True
maxX = 0
maxY = 0
ULpoints.append([maxX,maxY])
for patt in data:
    patternfull = patt.split(' ')
    if (len(patternfull) == 0): continue
    xy = patternfull[0] # (3,1)
    X,Y = map(int,xy[1:-1].split(','))
    if (X > maxX): maxX = X
    if (Y > maxY): maxY = Y
    if (Y == 0): inUL = False
    if inUL:
        ULpoints.append([X,Y])
    else:
        LRpoints.append([X,Y])
LRpoints.append([maxX,maxY])

# maxBla = [val,X,Y]
ULmaxChi2 = [0,0,0]
ULmaxInfGain = [0,0,0]
for (X,Y) in ULpoints:
    t = chi2(maxX,X,maxY,Y)
    if (t > ULmaxChi2[0]): ULmaxChi2 = [t,X,Y]
    t = infgain(maxX,X,maxY,Y)
    if (t > ULmaxInfGain[0]): ULmaxInfGain = [t,X,Y]
print "Upper Left hull (%s,%s), chi2=%s"%(ULmaxChi2[1],ULmaxChi2[2],ULmaxChi2[0])
print "Upper Left hull (%s,%s), InfGain=%s"%(ULmaxInfGain[1],ULmaxInfGain[2],ULmaxInfGain[0])

LRmaxChi2 = [0,0,0]
LRmaxInfGain = [0,0,0]
for (X,Y) in LRpoints:
    t = chi2(maxX,X,maxY,Y)
    if (t > LRmaxChi2[0]): LRmaxChi2 = [t,X,Y]
    t = infgain(maxX,X,maxY,Y)
    if (t > LRmaxInfGain[0]): LRmaxInfGain = [t,X,Y]
print "Lower Right hull (%s,%s), chi2=%s"%(LRmaxChi2[1],LRmaxChi2[2],LRmaxChi2[0])
print "Lower Right hull (%s,%s), InfGain=%s"%(LRmaxInfGain[1],LRmaxInfGain[2],LRmaxInfGain[0])


# A straightforward use of gnuplot.  The `debug=1' switch is used
# in these examples so that the commands that are sent to gnuplot
# are also output on stderr.
g = Gnuplot.Gnuplot(debug=0)
g.title('Convex Hull of %s (blue: Chi2, mangeta InfGain)'%getname(datafile))
g('set multiplot')
g('set noautoscale')
g('set xrange [0:%i]'%maxX)
g('set yrange [0:%i]'%maxY)
g('set arrow from %i,%i to %i,%i lt 3 lw 2'%(ULmaxChi2[1],ULmaxChi2[2]-(maxY/10),ULmaxChi2[1],ULmaxChi2[2]))
g('set arrow from %i,%i to %i,%i lt 4 lw 2'%(ULmaxInfGain[1],ULmaxInfGain[2]+(maxY/10),ULmaxInfGain[1],ULmaxInfGain[2]))
g('set arrow from %i,%i to %i,%i lt 3 lw 2'%(LRmaxChi2[1],LRmaxChi2[2]-(maxY/10),LRmaxChi2[1],LRmaxChi2[2]))
g('set arrow from %i,%i to %i,%i lt 4 lw 2'%(LRmaxInfGain[1],LRmaxInfGain[2]+(maxY/10),LRmaxInfGain[1],LRmaxInfGain[2]))
g('set style data linespoints') # give gnuplot an arbitrary command
g('set pointsize 5') # give gnuplot an arbitrary command
# Plot a list of (x, y) pairs (tuples or a numpy array would
# also be OK):
g.plot(ULpoints)
g.plot(LRpoints)

raw_input('Please press return to continue...\n')
