"""
myims module by David Floyd 31/07/2012
Useful tools for manipulating multiple (e.g. all-sky) images
Uses kd-trees for nearest neighbour search (in scipy.spatial)
You generate a tree using "gen_img_tree_wcs"
You then find which object is on which image using "match_img_wcs"
You can access a given image and return a flux at a given pixel using get_imfluxIf you wish to do this for many images, use get_imflux_loop
Some of these scripts require an input in the form of an ATPY table wiut hthe correct column headers (ra, dec).
"""

import numpy
from myutils import *
import pyfits, pywcs
import astrolib.coords as coords
import ds9
from scipy import *
from math import *
from scipy import spatial


def gen_img_tree_wcs(img_list,use_crval=True,proj_Euclid=False):
    map_ctrs = []
    for file in img_list:
        print ""
        print "%s:" % file
        hdulist = pyfits.open(file)
        if hdulist[0].header['CTYPE1'] == 'EQU--CAR':
            hdulist[0].header['CTYPE1'] = 'RA---CAR'
            hdulist[0].header['CTYPE2'] = 'DEC--CAR'
        wcs=pywcs.WCS(hdulist[0].header,naxis=2)
        if use_crval:
            map_ctr = wcs.wcs.crval
            print "Using CRVAL: ", map_ctr
        else:
            map_ctr = array([wcs.wcs_pix2sky(wcs.naxis2/2.,wcs.naxis1/2.,1)[0][0],wcs.wcs_pix2sky(wcs.naxis2/2.,wcs.naxis1/2.,1)[1][0]])
            print "Using calculated centre: ", map_ctr
        print "%f %f" % (map_ctr[0],map_ctr[1])
        equi = wcs.wcs.equinox
        if equi == 1950:
            print 'Equinox 1950: Precessing to J2000...'
            pos = coords.Position((map_ctr[0],map_ctr[1]),'b1950')
            map_ctr = pos.j2000()
        elif equi ==2000:
            print "J2000.0: no precession" 
        else:
            raise EquinoxError, "Don't know what to do with equi = %s : Not implemented yet!!!" % str(equi)
        if proj_Euclid:
            print "Euclidean projection -- multiplying RA by cos(dec)"
            # MULTIPLY RA by cos (dec) to get linear coord system!!!
            map_ctrs.append([map_ctr[0]*numpy.cos(numpy.radians(map_ctr[1])),map_ctr[1]])
        else:
            map_ctrs.append([map_ctr[0],map_ctr[1]])
        print map_ctrs[-1]
        #     if use_crval:
        #         map_ctrs.append(wcs.wcs.crval)
        #         #map_ctrs.append([wcs.wcs.crval[0]*numpy.cos(numpy.radians(wcs.wcs.crval[1])), wcs.wcs.crval[1]])
        #         #print [wcs.wcs.crval[0]*numpy.cos(numpy.radians(wcs.wcs.crval[1])), wcs.wcs.crval[1]]
        #     else:
        #         map_ctrs.append([map_ctr[0][0],map_ctr[1][0]])
        #         #map_ctrs.append([map_ctr[0][0]*numpy.cos(numpy.radians(map_ctr[1][0])),map_ctr[1][0]])
        #         #print [map_ctr[0][0]*numpy.cos(numpy.radians(map_ctr[1][0])),map_ctr[1][0]]
        hdulist.close()
    map_ctrs=numpy.array(map_ctrs)
    map_tree = spatial.KDTree(map_ctrs)
    print 'Returned map is in J2000'
    return map_tree,map_ctrs

def disp_img(hdulist):
    # DISPLAY:
    d = ds9.ds9()
    rslt=d.set_pyfits(hdulist)
    return d

def mark_img(d,p,text,regsz,color='green'):
    # DISPLAY:
    cmd_circ='regions command {circle '+str(p[0]+1)+' '+str(p[1]+1)+' '+str(regsz)+' # color='+color+'}'
    cmd_txt ='regions command {text ' +str(p[0]+1+regsz)+' '+str(p[1]+1+regsz+6)+' '+'# color='+color+' text="'+text+'" font="times 12 bold"}'
    rslt=d.set(cmd_circ)
    rslt=d.set(cmd_txt)
    print cmd_txt

def cutout(pt,filename,name,cutsz=256,disp_regsz=15,imsz=1024):
    """Usage... need to set these things up...
    #objno=[i for i,val in enumerate(twomass_rs.CATID) if val.startswith('23595879+0042066')]
    #pt=twomass_points[objno]
    #mapno=twomass_rs.NVSS_map_no[objno]
    #img=pyfits.open(nvss_img_list[mapno][0])
    #name=twomass_rs.CATID[objno][0]
    """
    hhh=pyfits.open(filename[0])
    d=myims.disp_img(hhh)
    wcs = pywcs.WCS(hhh[0].header,naxis=2)
    #scidata = hhh[0].data[0][0]
    scidata = numpy.ma.masked_array(hdulist[0].data[0][0],numpy.isnan(hdulist[0].data[0][0]))[:imsz,:imsz]
    pos = wcs.wcs_sky2pix(pt,0)[0] # data is 0-indexed! 
    cutout=scidata[int(pos[0]-cutsz/2):int(pos[0]+cutsz/2),int(pos[1]-cutsz/2):int(pos[1]+cutsz/2)]
    myims.mark_img(d,pos,name,disp_regsz)
    hhh.close()

def match_img_wcs(table, img_tree, tol,k=1,proj_Euclid=False,imsz=2.):
    print ""
    print "NOTE: match_img needs an ATpy table object with ra and dec columns"
    print ""
    print "Finding %i nearest neighbours..." % k
    print ""
    print "Assuming image half-side of %f degrees" % imsz
    print ""
    print ""
    if proj_Euclid:
        print "Using Euclidean projection.... multiplying RA by cos(dec)"
        points = numpy.array(zip(table.ra*numpy.cos(numpy.radians(table.dec)),table.dec))
    else:
        points = numpy.array(zip(table.ra,table.dec))
    print " Querying tree..."
    d_img,i_img = img_tree.query(points,k=k,distance_upper_bound=tol)
    # ####################################################################################################
    # You also need to check that any points near 360 degrees are caught by the appropriate tree point,...
    # ####################################################################################################
    print 'Checking tree coordinate wrap at 360 degrees, assuming image half-side %f degrees' % imsz
    if proj_Euclid:
        #edgeprox = [i for i, val in enumerate(points[:,0]) if val/cos(radians(points[i,1])) > edge]
        print "!!!!!!!Edge calculation not implemented for proj_Euclid...."
    else:
        print "Finding all objects withion %f degrees of 24h RA at all declinations...." % imsz
        edgeprox = [i for i, val in enumerate(table.ra) if val >= 360.-(imsz/(numpy.cos(numpy.radians(table.dec[i]))))]
    print 'Found %i points' % len(edgeprox)
    edge_points = points[edgeprox]
    # subtract 360 from these to make them have negative RA's and thus find a match to the 0h images!
    edge_points[:,0]-=360.
    d_img2,i_img2 = img_tree.query(edge_points,k=k,distance_upper_bound=tol)
    d_img[edgeprox] = d_img2
    i_img[edgeprox] = i_img2
    # ####################################################################################################
    print "Still need to correct hi latitude objects.."
    return d_img, i_img


def get_imflux_loop(img_list, table, points, IDs, i_map, tol, pause=False, display=True, use_max=False, amelia=False, use_tot=False, verbose=False, imsz=1024, disp_regsz=5, sub_frame='', div_frame='', use_square=False,use_square_bg=False):
    # tol is PIXEL RANGE to search over OR (if use_tot = True) then the RADIUS of circular aperture (in PIXELS) to sum over...
    if use_tot:
        print "Setting up circular aperture of radius %i pixels, centred at centre of image" % tol
        # Define a basic circle for translation to any position we like...
        circ_base = numpy.array([(x,y) for x in arange(imsz) for y in arange(imsz) if (x-imsz/2+1)**2+(y-imsz/2+1)**2 <= tol**2])
    tile_flux = numpy.array(zeros(len(points)))
    tile_mean = numpy.array(zeros(len(points)))
    tile_stdv = numpy.array(zeros(len(points)))
    # Index map number with m:::
    for m, mapn in enumerate(img_list):
        #m=91
        #mapn=nvss_img_list[m]
        if verbose: print mapn
        hdulist = pyfits.open(mapn)
        if hdulist[0].header['CTYPE1'] == 'EQU--CAR':
            hdulist[0].header['CTYPE1'] = 'RA---CAR'
            hdulist[0].header['CTYPE2'] = 'DEC--CAR'
        wcs = pywcs.WCS(hdulist[0].header,naxis=2)
        naxis = hdulist[0].header['NAXIS']
        # Define list of points that are in this image...
        selection = [i for i, mapno in enumerate(i_map) if mapno == m]
        points_to_query = points[selection]
        names = IDs[selection]
        if verbose: print "Selection: ", names, points_to_query
        no_in_map = len(selection)
        # DISPLAY?
        if display or verbose: d = disp_img(hdulist)
        if no_in_map > 0:
            # Set up list of points... and precess if necessary...
            new_points_to_query = points_to_query - points_to_query
            equi = wcs.wcs.equinox
            if equi == 1950: # need to precess points back to 1950 to match WCS!
                print ' '
                print 'Uh-oh, it looks like your WCS is in B1950!'
                print 'Assuming input points are J2000...'
                print 'Now precessing your input points to B1950 to match WCS...'
                for p, point in enumerate(points_to_query):
                    pos  = coords.Position(point)
                    newpos = pos.b1950()
                    new_points_to_query[p][0] = newpos[0]
                    new_points_to_query[p][1] = newpos[1]
                print 'PRECESSED points to B1950 to match WCS!!!'
                print ' '
            else:
                print ' '
                print 'WCS is in J2000'
                print 'Assuming input points are also J2000...'
                print 'No precession!'
                print ' '
                new_points_to_query = points_to_query
            print 'Map %i: %s    No. of objects: %i' % (m, mapn, no_in_map)
            # Now read in the science data. How you do this depends on NAXIS
            if naxis == 2:
                #scidata = hdulist[0].data
                scidata = numpy.ma.masked_array(hdulist[0].data,numpy.isnan(hdulist[0].data))[:imsz,:imsz]
                # THIS IS SPECIAL FOR ROSAT RASS IMAGES: need to divide thru by exptime (and possibly sutract background!)
                if len(sub_frame) == len(img_list):
                    hdusub = pyfits.open(sub_frame[m])
                    subdata=numpy.ma.masked_array(hdusub[0].data,numpy.isnan(hdusub[0].data))[:imsz,:imsz]
                    scidata-=subdata
                    hdusub.close()
                if len(div_frame) == len(img_list):
                    print 'Divinding by %s' % div_frame[m]
                    hdudiv = pyfits.open(div_frame[m])
                    divdata=numpy.ma.masked_array(hdudiv[0].data,numpy.isnan(hdudiv[0].data))[:imsz,:imsz]
                    scidata/=divdata
                    hdudiv.close()
            elif naxis == 3:
                #scidata = hdulist[0].data[0]
                scidata = numpy.ma.masked_array(hdulist[0].data[0],numpy.isnan(hdulist[0].data[0]))[:imsz,:imsz]
            elif naxis == 4:
                #scidata = hdulist[0].data[0][0]
                scidata = numpy.ma.masked_array(hdulist[0].data[0][0],numpy.isnan(hdulist[0].data[0][0]))[:imsz,:imsz]
            map_mean = numpy.mean(scidata)
            map_std = numpy.std(scidata)
            if map_mean >= 0.01 or map_std >= 0.01:
                print "NOISY MAP: %i, %s. Mean level %f +/- %f. Treat with caution!!!" % (m, mapn, map_mean, map_std)
                col = "blue"
            else:
                col = "green"
            pos1 = wcs.wcs_sky2pix(new_points_to_query,0) # data is 0-indexed! 
            for j, p1 in enumerate(pos1):
                pixpos = [iround(p1[1]),iround(p1[0])] # AND x-y inverted!!
                if verbose:
                    print "sky position  : ",new_points_to_query[j]
                    print "pixel position: ",p1
                    print pixpos
                #if 0<pixpos[0]<wcs.naxis2 and 0<pixpos[1]<wcs.naxis1:
                if 0<=pixpos[0]<=imsz-1 and 0<=pixpos[1]<=imsz-1:
                    if amelia:
                        # Do Amelia's flux-gettiung....
                        print "Amelias flux-getting stuff not yet implemented"
                        tile_flux[selection[j]] = 0.0 
                        #tile_flux[selection[j]] = numpy.sum(scidata[pixpos[0]-tol/2:pixpos[0]+tol/2,pixpos[1]-tol/2,pixpos[1]+tol/2])
                    elif use_max:
                        # Find MAXIMUM pixel value within tol pixel of your position
                        tile_flux[selection[j]] = getmaxval(scidata,pixpos,tol)
                    elif use_tot:
                        #tile_flux[selection[j]] = gettotval(scidata,pixpos,tol)
                        #circ = numpy.array([(x,y) for x in arange(imsz) for y in arange(imsz) if (x-pixpos[0])**2+(y-pixpos[1])**2 <= rad**2])
                        circ = circ_base.copy()
                        circ[:,0]-=(imsz/2-pixpos[0])
                        circ[:,1]-=(imsz/2-pixpos[1])
                        try:
                            tile_flux[selection[j]] = numpy.sum(scidata[[circ[:,0]],[circ[:,1]]])
                        except IndexError, e:
                            print "IndexError", e, " in gettotval"
                            print "Object %i too close to edge of map %i: %s" % (j,m,mapn)
                            print "Setting BAD"
                            tile_flux[selection[j]] = -999.
                    elif use_square_bg:
                        thick=2
                        if verbose:
                            print "Applying square apertures... with background subtraction"
                        data_square1=square_ap(pixpos[0],pixpos[1],tol)
                        data_square2=square_ap(pixpos[0],pixpos[1],tol+2*thick)
                        try:
                            tf1=numpy.sum(scidata[data_square1])
                            tf2=numpy.sum(scidata[data_square2])
                            ann=(tf2-tf1)/(len(data_square2)-len(data_square1))
                            box=tf1/len(data_square1)
                            tile_flux[selection[j]] = box-ann
                        except IndexError, e:
                            print "IndexError", e, " in square annulus calc..."
                            print "Object %i too close to edge of map %i: %s" % (j,m,mapn)
                            print "Setting BAD"
                            tile_flux[selection[j]] = -999.
                    else:
                        # Just get the flux at this pixel position...
                        tile_flux[selection[j]] = getpixval(scidata,pixpos)
                    if display or verbose:
                        mark_img(d,p1,names[j],disp_regsz,color=col)
                else:
                    print "OBJECT %i OFF EDGE OF TILE %s !!!" % (selection[j], img_list[m])
                    tile_flux[selection[j]] = -999.
                tile_mean[selection[j]] = map_mean
                tile_stdv[selection[j]] = map_std
        else:
            print 'No objects in map %s' % (img_list[m])
        if pause or verbose:
            raw_input("push enter to continue")
        hdulist.close()
    return tile_flux, tile_mean, tile_stdv


def get_imflux(mapn, points, IDs, display=True,disp_regsz = 5):
    # disp_regsz =  radius of circle to display over each source...
    tol = 0 # pixel range to loop for peak pixel over...
    tile_flux = numpy.array(zeros(len(points)))
    print "opening %s" % mapn
    hdulist = pyfits.open(mapn)
    if display: 
        d = disp_img(hdulist)
    print "Obtaining WCS and data..."
    wcs = pywcs.WCS(hdulist[0].header,naxis=2)
    naxis = hdulist[0].header['NAXIS']
    if naxis == 2:
        scidata = hdulist[0].data
    elif naxis == 3:
        scidata = hdulist[0].data[0]
    elif naxis == 4:
        scidata = hdulist[0].data[0][0]
    else:
        print "Dunno how to deal with FITS file with NAXIS = %i" % naxis
    print "Data has shape: ", shape(scidata)
    #scidata = numpy.ma.masked_array(hdulist[0].data[0][0],numpy.isnan(hdulist[0].data[0][0]))
    #mean = numpy.mean(scidata)
    #stdev = numpy.std(scidata)
    pos = wcs.wcs_sky2pix(points,0) # data is 0-indexed! 
    for j, p in enumerate(pos): 
        pixpos = [iround(p[1]),iround(p[0])] # AND x-y inverted!!
        if 0<pixpos[0]<wcs.naxis2-1 and 0<pixpos[1]<wcs.naxis1-1:
            flux = getpixval(scidata,pixpos)
            #flux_max = getmaxval(scidata,pixpos,tol)
            print '    ', j, IDs[j], pixpos, flux #,flux_max
            tile_flux[j] = flux
        else:
            print "OBJECT %i: %s OFF EDGE OF TILE!!!" % (j, IDs[j])
            tile_flux[j] = nan
        print "Tile flux: %f" % tile_flux[j]
        if display:
            print "Marking object %i " % j
            if tile_flux[j] > 5.:
                mark_img(d,p,IDs[j],disp_regsz,color="red")
            elif tile_flux[j] > 0.6:
                mark_img(d,p,IDs[j],disp_regsz,color="blue")
            else:
                 mark_img(d,p,IDs[j],disp_regsz)
    hdulist.close()
    return tile_flux#,mean,stdev



