"""
myfuns module by David Floyd 31/07/2012
Useful miscellaneous functions
Some of these scripts require an input in the form of an ATPY table with hthe correct column headers (e.g. ra, dec, CATID).

"""

import numpy, pyfits, sys, os
from scipy import spatial, nan
import pkg_resources

# Directory walk, using generators...
def dirwalk(dir):
    for f in os.listdir(dir):
        fullpath = os.path.join(dir,f)
        if os.path.isdir(fullpath) and not os.path.islink(fullpath):
            for x in dirwalk(fullpath):  # recurse into subdir
                yield x
        else:
            yield fullpath

def load_just_installed(name):
# This neat trick will load  something newly installed
    pkg_resources.get_distribution(name).activate()
    print "Now you can import %s" % name
    return

def iround(x):
    """iround(number) -> integer
    Round a number to the nearest integer."""
    return int(round(x) - .5) + (x > 0)

def getpixval (img,point):
    # return value at given pixel position
    try:
        val = img [point[0],point[1]]
    except ValueError, e:
        print "ValueError", e, " in getpixval"
        print "Setting pixel value at ", point, " to nan"
        val = nan
    return val

def getmaxval(data,pixpos,tol):
    # Return maximum pixvalue within some region (defined by tol) around pixpos:
    try:
        val = (data[pixpos[0]-tol:pixpos[0]+tol,pixpos[1]-tol:pixpos[1]+tol]).max()
    except ValueError, e:
        print "ValueError", e, " in getmaxval"
        print pixpos[0]-tol, pixpos[0]+tol,pixpos[1]-tol, pixpos[1]+tol
        print "TRYING GET_PIXVAL..."
        val = getpixval(data,pixpos)
        print "THAT WORKED!"
    return val

def gettotval(data,pixpos,rad,imsz=1024):
    circ=numpy.array([(x,y) for x in arange(imsz) for y in arange(imsz) if (x-pixpos[0])**2+(y-pixpos[1])**2 <= rad**2])
    val = numpy.sum(scidata[[circ[:,0]],[circ[:,1]]])
    return val

def do_ap_phot(img,point,region):
    # Even better: do some sort of aperture photom?...TBA
    print "do_ap_phot NOT YET IMPLEMENTED!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
    val = 0.
    return val

#def findmap(points1,points2):
#    for point2 in points2:
#        for point1 in points1:
#            dist=coord.ang_sep(point1[0],point1[1],point2[0],point2[1])
#            print dist

def match_string_lists(list1, list2):
    set1=set(list1)
#    matches = [i for i, item in enumerate(list1) if item in set2]
    for obj in set1:
        match = index(list1,)
#[i for i, item in enumerate(list1) if item in set2]
    return matches

# insideabox: returns boolean array with an element for each point (x,y) that falls inside the box ([x,y],[x,y],[x,y],[x,y])
def insideabox(points,box):
    if box[2,0] > box [0,0]:
        box[2,0]-=360.
    inside = []
    for point in points:
        xin = box[2,0]+.1 <= point[0] <= box[0,0]-.1
        yin = box[0,1]+.1 <= point[1] <= box[2,1]-.1
        #        print point, xin, yin
        if xin and yin:
            inside.append(True)
        else:
            inside.append(False)
    return numpy.array(inside)

# BENCHMARK RUN::: TBA???
# USE timeit instead:
# http://docs.python.org/library/timeit.html
#def run_bm(task):
#    st=time.time()
#    task
#    et=time.time()
#    print "time taken for %s = %f seconds" % ( task, et-st )

def match_coords(coords1, coords2, tol):
# Create tree from coords1 
    tree = spatial.KDTree(coords1)
# Lookup neighbours to each point, by point in coords2:
    q = tree.query_ball_point(coords2,tol)
#    d_map,q = tree.query(coords2,k=1,distance_upper_bound=tol)
    index1=[]
    for qq in q:
        try:
            ans=qq[0]
        except IndexError, e:
            ans=-1
        index1.append(ans)
    index1=numpy.array(index1)
    keeprows2 = numpy.array([i for i, ind in enumerate(index1) if ind!=-1])
    return index1,keeprows2

def index(vec, statement):
    return [i for i,el in enumerate(vec) if statement]

def findname(table,names):
    print """NOTE: findname requires an input atpy table object with a column called CATID (containing the names of the sources)"""
    index=[]
    for j,name in enumerate(names):
        try:
            index.append([i for i,ID in enumerate(table.CATID) if ID[:len(name)] == name][0])
        except IndexError:
            print "# %i: object %s not found " % (j, name)
    return index

# dec2deg, ra2deg Written by Enno Middelberg 2001
def dec2deg(s):
    print s
    inp = s.split(":")
    print inp
    if len(inp)<3:
        raise ValueError("Too few fields in H:M:S string.")
    elif len(inp) > 3:
        raise ValueError("Too many fields in H:M:S string.")
    hh=abs(float(inp[0]))
    mm=float(inp[1])/60.
    ss=float(inp[2])/3600.
    if float(inp[0]) < 0 or inp[0]=='-00' or inp[0]=='-0':
        print" That's -"+str(hh+mm+ss)+" degrees"
        dec = -(hh+mm+ss)
    else:
        print" That's "+str(hh+mm+ss)+" degrees"
        dec = hh+mm+ss
    return (dec)

def ra2deg(s):
    print s
    inp = s.split(":")
    if len(inp)<3:
        raise ValueError("Too few fields in H:M:S string.")
    elif len(inp) > 3:
        raise ValueError("Too many fields in H:M:S string.")
    hh=(float(inp[0]))*15.
    mm=(float(inp[1])/60.)*15.
    ss=(float(inp[2])/3600.)*15.
    print" That's "+str(hh+mm+ss)+" degrees"
    ra = hh+mm+ss
    return (ra)

def str_HMS(s):
    """Convert a string "H:M:S" to a tuple of numbers (H, M, S).
    H and M must be integers, S may be a float.
    """
    L = s.split(":")
    if len(L) < 3:
        raise ValueError("Too few fields in H:M:S string.")
    elif len(L) > 3:
        raise ValueError("Too many fields in H:M:S string.")
    H = int(L[0])
    M = int(L[1])
    S = float(L[2])
    return (H, M, S)

def HMS_Seconds(hours, minutes, seconds):
    """Convert hours minutes seconds to seconds."""
    return hours*60*60 + minutes*60 + seconds

def HMS_Hours(hours, minutes, seconds):
    """Convert hours minutes seconds to hours."""
    return hours + minutes/60.0 + seconds/(60.0**2)

def Hours_HMS(t):
    """Convert time t in hours to hours minutes seconds."""
    hours = int(t)
    t = (t - hours)*60
    minutes = int(t)
    seconds = (t - minutes)*60
    return (hours, minutes, seconds)

def Seconds_HMS(t):
    """Convert time t in seconds to hours minutes seconds."""
    hours, t = divmod(t, 60*60)
    minutes, seconds = divmod(t, 60)
    return (hours, minutes, seconds)
