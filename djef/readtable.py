#filelist=glob.glob('/Volumes/Work/floyd/DB/2MASS/psc_*')
#filelist=glob.glob('/Volumes/Work/floyd/DB/2MASS/xsc_*')

#tset=atpy.TableSet()

#    for file in filelist:
#    thisfile=open( file ,'r')
#    theselines=thisfile.readlines()
#    thisfile.close()
#    lines=[lines,theselines]
# The above method constructs lines quickly, but must then process them...
# ... or use atpy: for ascii this just calls asciitable...
# this is slower but creates a structure.
#        tset.append(t)
# Define a Class for storing data....

import atpy,numpy
def readtable(file,type):
    class cat:
        pass
    t=atpy.Table(file,type=type)
    col_list = t.keys()
    for col in col_list:
        vars( cat )[ col ] = []
    for i, col in enumerate( col_list ):
        for entry in t: 
            vars( cat )[ col ].append( entry[ i ] )
    for col in col_list:
        vars( cat )[ col ] = numpy.array( vars( cat )[ col ] )
    return cat


