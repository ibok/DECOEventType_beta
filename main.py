#!/usr/bin/env python

"""Runs through a given folder of images and prints out the type of event for
   each blob group in that event.
   Stats: Weighted accuracy is roughly 85%, Processing rate is 2 sec/image

   For ease of parsing I/O, Output is in the format: id:event_type.

   Appends first letter of suspected event type to name of file."""

import argparse, math, os, sys
import numpy as np
from PIL import Image
from skimage import measure

class Line:
    """Line class that builds a line 'object' that is able to run various
       mathematical tests on itself"""

    def __init__(self, m, p):
        self.m = float(m)
        self.p = (float(p[0]), float(p[1]))

    def calculate_y(self, x):
        return self.m*(x-self.p[0]) + self.p[1]

    def intersect(self, line):
        if line.m == self.m:
            return NaN
        else:
            xpos = (line.p[1] - self.p[1] + self.m*self.p[0] + line.m*line.p[0])/(self.m - line.m)
            return (xpos, self.calculate_y(xpos))

class Blob:
    """Class that defines a 'blob' in an image: the contour of a set of pixels
       with values above a given threshold."""

    def __init__(self, x, y):
        """Define a counter by its contour lines (an list of points in the xy
           plane), the contour centroid, and its enclosed area."""
        self.x = x
        self.y = y
        self.xc = np.mean(x)
        self.yc = np.mean(y)

        # Find the area inside the contour
        self.area = 0.
        n = len(x)
        for i in range(0, n):
            self.area += 0.5*(y[i]+y[i-1])*(x[i]-x[i-1])

        # Calculate perimeter of blob using cycling
        self.perimeter = 0
        end = len(self.x) - 1
        for i in range(0, end):
            a = (self.x[i], self.y[i])
            b = (self.y[i+1], self.y[i+1])
            self.perimeter += distance(a, b)

        # Find furthest distance between any two blob contour points
        self.maxdist = 0
        end = len(self.x)
        for i in range(0, end):
            for j in range(i, end):
                pt1 = (self.x[i], self.y[i])
                pt2 = (self.y[j], self.y[j])
                pdist = distance(pt1, pt2)
                if pdist > self.maxdist:
                    self.maxdist = pdist

    def distance(self, blob):
        """Calculate the distance between the centroid of this blob contour and
           another one in the xy plane."""
        return np.sqrt((self.xc - blob.xc)**2 + (self.yc-blob.yc)**2)

def distance(a, b):
    return math.sqrt((a[1] - b[1])**2 + (a[0] - b[0])**2)

class BlobGroup:
    """A list of blobs that is grouped or associated in some way, i.e., if
       their contour centroids are relatively close together."""

    def __init__(self):
        """Initialize a list of stored blobs and the bounding rectangle which
        defines the group."""

        self.blobs = []
        self.count = 0
        self.area = 0.
        self.b_area = 0.
        self.perimeter = 0.

        self.xmin =  1e10
        self.xmax = -1e10
        self.ymin =  1e10
        self.ymax = -1e10

    def addBlob(self, blob):
        """Add a blob to the group and enlarge the bounding rectangle of the
           group."""
        self.blobs.append(blob)
        self.xmin = min(self.xmin, blob.x.min())
        self.xmax = max(self.xmax, blob.x.max())
        self.ymin = min(self.ymin, blob.y.min())
        self.ymax = max(self.ymax, blob.y.max())
        self.cov  = None
        self.count += 1
        self.b_area += blob.area
        self.perimeter += blob.perimeter

    def getBoundingBox(self):
        """Get the bounding rectangle of the group."""
        return (self.xmin, self.xmax, self.ymin, self.ymax)

    def getSquareBoundingBox(self):
        """Get the bounding rectangle, redefined to give it a square aspect
           ratio."""
        xmin, xmax, ymin, ymax = (self.xmin, self.xmax, self.ymin, self.ymax)
        xL = np.abs(xmax - xmin)
        yL = np.abs(ymax - ymin)
        if xL > yL:
            ymin -= 0.5*(xL-yL)
            ymax += 0.5*(xL-yL)
        else:
            xmin -= 0.5*(yL-xL)
            xmax += 0.5*(yL-xL)
        return (xmin, xmax, ymin, ymax)

    def getSubImage(self, image):
        """Given an image, extract the section of the image corresponding to
           the bounding box of the blob group."""
        ny,nx = image.shape
        x0,x1,y0,y1 = bg.getBoundingBox()
        self.area = (x1 - x0)*(y1 - y0)

        # Account for all the weird row/column magic in the image table...
        i0,i1 = [ny - int(t) for t in (y1,y0)]
        j0,j1 = [int(t) for t in (x0,x1)]

        # Add a pixel buffer around the bounds, and check the ranges
        buf = 1
        i0 = 0 if i0-buf < 0 else i0-buf
        i1 = ny-1 if i1 > ny-1 else i1+buf
        j0 = 0 if j0-buf < 0 else j0-buf
        j1 = nx-1 if j1 > nx-1 else j1+buf

        return image[i0:i1, j0:j1]

    def getRawMoment(self, image, p, q):
        """Calculate the image moment given by
           M_{ij}=\sum_x\sum_y x^p y^q I(x,y)
           where I(x,y) is the image intensity at location x,y."""
        nx,ny = image.shape
        Mpq = 0.
        if p == 0 and q == 0:
            Mpq = np.sum(image)
        else:
            for i in range(0,nx):
                x = 0.5 + i
                for j in range(0,ny):
                    y = 0.5 + j
                    Mpq += x**p * y**q * image[i,j]
        return Mpq

    def getCovariance(self, image):
        """Get the raw moments of the image region inside the bounding box
           defined by this blob group and calculate the image covariance
           matrix."""
        if self.cov == None:
            subImage = self.getSubImage(image).transpose()
            M00 = self.getRawMoment(subImage, 0, 0)
            M10 = self.getRawMoment(subImage, 1, 0)
            M01 = self.getRawMoment(subImage, 0, 1)
            M11 = self.getRawMoment(subImage, 1, 1)
            M20 = self.getRawMoment(subImage, 2, 0)
            M02 = self.getRawMoment(subImage, 0, 2)
            xbar = M10/M00
            ybar = M01/M00
            self.cov = np.vstack([[M20/M00 - xbar*xbar, M11/M00 - xbar*ybar],
                                  [M11/M00 - xbar*ybar, M02/M00 - ybar*ybar]])
        return self.cov

    def getPrincipalMoments(self, image):
        """Return the maximum and minimum eigenvalues of the covariance matrix,
           as well as the angle theta between the maximum eigenvector and the
           x-axis."""
        cov = self.getCovariance(image)
        u20 = cov[0,0]
        u11 = cov[0,1]
        u02 = cov[1,1]

        theta = 0.5 * np.arctan2(2*u11, u20-u02)
        l1 = 0.5*(u20+u02) + 0.5*np.sqrt(4*u11**2 + (u20-u02)**2)
        l2 = 0.5*(u20+u02) - 0.5*np.sqrt(4*u11**2 + (u20-u02)**2)
        return l1, l2, l1/l2, theta

def findBlobs(image, threshold, minArea=2.):
    """Pass through an image and find a set of blobs/contours above a set
       threshold value.  The minArea parameter is used to exclude blobs with an
       area below this value."""
    blobs = []
    ny, nx = image.shape

    # Find contours using the Marching Squares algorithm in the scikit package.
    contours = measure.find_contours(image, threshold)
    for contour in contours:
        x = contour[:,1]
        y = ny - contour[:,0]
        blob = Blob(x, y)
        if blob.area >= minArea:
            blobs.append(blob)
    return blobs

def groupBlobs(blobs, maxDist):
    """Given a list of blobs, group them by distance between the centroids of
       any two blobs.  If the centroids are more distant than maxDist, create a
       new blob group."""
    n = len(blobs)
    groups = []
    if n >= 1:
        # Single-pass clustering algorithm: make the first blob the nucleus of
        # a blob group.  Then loop through each blob and add either add it to
        # this group (depending on the distance measure) or make it the
        # nucleus of a new blob group
        bg = BlobGroup()
        bg.addBlob(blobs[0])
        groups.append(bg)

        for i in range(1, n):
            bi = blobs[i]
            isGrouped = False
            for group in groups:
                # Calculate distance measure for a blob and a blob group:
                # blob just has to be < maxDist from any other blob in the group
                for bj in group.blobs:
                    if bi.distance(bj) < maxDist:
                        group.addBlob(bi)
                        isGrouped = True
                        break
            if not isGrouped:
                bg = BlobGroup()
                bg.addBlob(bi)
                groups.append(bg)

    return groups

def get_info(group):
    # Returns parsed image data of any picture
    tearr = []
    maxdist = [0]
    tarr = [0, 0, 0]
    intervs = 100.
    first = True

    for b in group.blobs:
        for i in range(0, len(b.x)):
            tearr.append([float(b.x[i]), float(b.y[i])])
    for i in range(0, len(tearr)):
        for j in range(i, len(tearr)):
            if (distance(tearr[i], tearr[j]) > maxdist[0]):
                maxdist = [distance(tearr[i], tearr[j]), tearr[i], tearr[j]]
    if maxdist[1][0] < maxdist[2][0]:
        leftmost = maxdist[1][0]
        y1 = maxdist[1][1]
        rightmost = maxdist[2][0]
    else:
        leftmost = maxdist[2][0]
        y1 = maxdist[2][1]
        rightmost = maxdist[1][0]
        first = False
    num = maxdist[1][1] - maxdist[2][1]
    dem = maxdist[1][0] - maxdist[2][0]
    if first:
        return [map(int, maxdist[1]), map(int, maxdist[2])]
    else:
        return [map(int, maxdist[2]), map(int, maxdist[1])]

def opt_contour(group):
    # Removes duplicate contour points
    tearr = []
    for b in group.blobs:
        for i in range(0, len(b.x)):
            tearr.append((float(b.x[i]), float(b.y[i])))
    return sorted(set(tearr))

def calcslp(num, dem):
    # Automatic slope filtering function that filters out inf or -inf slopes
    if dem == 0 and num == 0:
        return NaN
    elif dem == 0:
        return np.sign(num)*np.sqrt(sys.maxint)
    else:
        return float(num)/dem

def findDist(line, point):
    return abs(-line.m*point[0]+point[1]-(line.p[1] - line.m*line.p[0]))/float(np.sqrt((-line.m)**2 + 1))

def print_type(iid, type):
    print >>f, str(iid) + ',' + get_chars(type, 0)

def exit(): # stops program
    raise SystemExit

# Get an image file name from the command line
p = argparse.ArgumentParser(description="Histogram luminance of a JPG image")
p.add_argument("folder", nargs=1,
               help="JPG folder name")
p.add_argument("-x", "--threshold", dest="thresh", type=float,
               default=None,
               help="Threshold the pixel map")
blobGroup = p.add_argument_group("Blob Identification")
blobGroup.add_argument("-c", "--contours", dest="contours", type=float,
                       default=None,
                       help="Identify blob contours above some threshold")
blobGroup.add_argument("-d", "--distance", dest="distance", type=float,
                       default=75.,
                       help="Group blobs within some distance of each other")
blobGroup.add_argument("-a", "--min-area", dest="area", type=float,
                       default=2.,
                       help="Remove blobs below some minimum area")
args = p.parse_args()

# Recursively list all files in a directory, checking for valid image types
filegen = os.walk(args.folder[0])

flist = []
img_types = ('jpg jpeg bmp png eps gif im j2k j2p jpx msp pcx png ppm pgm pbm' +
            'spi tiff webp xbm xv').split(' ') #Most of these are useless...

# Create file list!
for root, afile, i in filegen:
    for fil in i:
        if fil.split('.')[-1] in img_types:
            flist.append((root, root + '/' + fil))

info_arr = []
pt_arr = []

# Check if the folder contains valid image files and if not, exit
if len(flist) == 0:
    print('No image files were found.') and exit()

# Now that we know the folder has images, write a classification file.
f = open('classifications.out', 'w')

def get_chars(x, opt): #Returns string version of type and/or appending char.
    if opt == 0:
        return {
            0 : 'null',
            1 : 'spot',
            2 : 'worm',
            3 : 'track',
            4 : 'ambig',
            5 : 'big_spot',
            6 : 'track_lowconf'
        }[x]
    else:
        return {
            '0' : '_x',
            '1' : '_s',
            '2' : '_w',
            '3' : '_t',
            '4' : '_a',
            '5' : '_b',
            '6' : '_l'
        }[x]

for files in flist:
    # Fix filename if it contains extra slash
    efile = files[1].replace('//','/')

    # Prepare image for filename change
    fullname = efile.split('/')[-1]
    iid = fullname.split('.')[0]
    tail = fullname.split('.')[-1]

    # Load each image and convert pixel values to grayscale intensities
    img = Image.open(efile).convert("L")
    image = []
    pix = img.load()

    # Stuff image values into a 2D table called "image"
    nx = img.size[0]
    ny = img.size[1]
    x0, y0, x1, y1 = (0, 0, nx, ny)
    for y in xrange(ny):
        image.append([pix[x, y] for x in xrange(nx)])

    # Find average pixel value (grayscale) for noise classification
    image = np.array(image, dtype=float)
    pixavg = sum(sum(image))/(len(image)*len(image[0]))
    if pixavg > 4: # Skips picture if av. px. val is too high
        type = 7
        print_type(iid, type)
        continue

    # Calculate contours using the scikit-image marching squares algorithm,
    # store as Blobs, and group the Blobs into associated clusters

    blobs = findBlobs(image, threshold=args.contours, minArea=args.area)
    groups = groupBlobs(blobs, maxDist=args.distance)

    # Go through each group
    for g in groups:
        info_arr.append(get_info(g))
        pt_arr.append(opt_contour(g))

    # Apply a threshold to the image pixel values
    if args.thresh != None:
        image[image < args.thresh] = 0.

    if args.contours != None:
        for i, bg in enumerate(groups):
            X0, X1, Y0, Y1 = bg.getSquareBoundingBox()
            l1, l2, r, theta = bg.getPrincipalMoments(image)
            ecc = (np.sqrt( l1**2 - l2**2 ) / l1)
            # Calculate summative distance from the maximum distance
            # line of all points in the group's blob's contours, then weight.
            init = info_arr[i][0]
            final = info_arr[i][1]
            _i = Line(calcslp((float(info_arr[i][0][1]) - info_arr[i][1][1]),(float(info_arr[i][0][0]) - info_arr[i][1][0])), init)
            tdist = 0.
            mlength = distance(init, final)
            for pt in pt_arr[i]:
                tdist += findDist(_i, pt)
            factr = tdist/mlength #arbitrary normalizing factor

            """Sequence of statements that calculates the weighted
            center of image contours (by area), then uses the centerpoint
            ratios in both the X and Y direction to find the centerpoint
            displacement product (cdrp) and difference (cdrd) of the image
            group."""

            bgtarea = 0.
            carr = [0, 0]
            for b in bg.blobs:
                bgtarea += b.perimeter
                carr[0] += b.xc*b.perimeter
                carr[1] += b.yc*b.perimeter
            centx = carr[0]/float(bg.count*bgtarea)
            centy = carr[1]/float(bg.count*bgtarea)
            diffx_a = centx - bg.xmin
            diffx_b = bg.xmax - centx
            diffy_a = centy - bg.ymin
            diffy_b = bg.ymax - centy
            ratx = min([diffx_a, diffx_b])/max([diffx_a, diffx_b])
            raty = min([diffy_a, diffy_b])/max([diffy_a, diffy_b])
            cdrp = ratx*raty
            cdrd = abs(ratx-raty)

            """Code block that analyzes the most likely type of event inside
            the image group currently being analyzed."""

            # 0 == null/noise ; 1 == Spot ; 2 == Worm ; 3 == Track ; 4 == Ambiguous;
            # 5 = Alpha Particle; 6 = Track, low confidence
            if args.contours == 40:
                type = 0
                if ecc > 0.99993 and l1 > 700:
                    type = 4
                elif bg.b_area < 4 or mlength < 6 or mlength < 13 and (r >= 0.2 and ecc < 0.945 or bg.b_area < 7) or ecc < 0.7:
                    if bg.b_area > 50:
                        if bg.b_area > 62 and mlength > 10.:
                            type = 2
                        else:
                            type = 5
                    else:
                        type = 1
                else:
                    if cdrd > 0.55:
                        type = 2
                    elif bg.b_area > 100 and (l1 > 100 and l1/10 > l2) and mlength > 30:
                        if factr > 9 or l1/5 > mlength and factr > 3.9 or mlength > 40 and mlength < 80 and bg.b_area > 100 and factr > 5:
                            type = 2
                        else:
                            type = 3
                    elif ecc > 0.9995 and mlength > 40 and cdrp > 0.8:
                        if ecc > 0.9998 and mlength > 90 and mlength < 130:
                            type = 3
                        else:
                            type = 4
                    else:
                        if (cdrp > 0.978 and cdrd < 0.01 or cdrp > 0.96 and cdrd < 0.0047 or cdrp > 0.9 and r < 0.02 and cdrd < 0.055 and factr < 5.) and ecc > 0.96:
                            if bg.b_area > 33:
                                type = 3
                            elif ecc > 0.999:
                                type = 2
                            else:
                                type = 2
                        elif ecc < 0.985 and r < 0.21 and (cdrd < 0.015 or cdrp > 0.88) and mlength > 9 and mlength < 18:
                            type = 3
                        elif ecc < 0.985 and ecc > 0.97 and r < 0.21 and mlength > 9 and mlength < 18 and cdrp > 0.83 and bg.b_area < 30:
                            type = 3
                        elif ecc > 0.975 and ecc < 0.99 and r < 0.22 and factr < 3.7 and mlength > 7.6 and mlength < 18 and cdrd < 0.023:
                            type = 3
                        elif ecc > 0.99 and l1 < 15. and cdrp > 0.86 and cdrd < 0.1 and bg.b_area > 28 and bg.b_area < 35:
                            type = 3
                        else:
                            if factr > 4.6 and bg.b_area < 100 or bg.b_area < 24 or ecc <= 0.978 or r > 0.2 and factr > 4.:
                                type = 2
                            else:
                                if l1 > 100 and l2 < 12 and bg.b_area > 40 and bg.b_area < 60 and factr < 3.8:
                                    type = 3
                                elif (abs(ratx) > 0.99 or abs(raty) > 0.99) and abs(ratx) > 0.93 and abs(raty) > 0.93 and ecc > 0.99 and factr < 2.95 and factr > 2. and cdrd > 0.05:
                                    type = 3
                                elif cdrd < 0.7 and factr > 6:
                                    type = 2
                                elif (cdrp > 0.9 and cdrd < 0.02 and factr < 3.1 and ecc > 0.993 and mlength > 12) and bg.b_area < 82: #random magic
                                    type = 3
                                elif ((cdrp > 0.6 and ecc > 0.9923 and factr < 3.1 or cdrp > 0.88) and (cdrd < 0.03 or abs(ratx) > 0.996 or abs(raty) > 0.996) and bg.b_area < 100):
                                    type = 3
                                else:
                                    if ecc > 0.999 and (l1 > 90 and l2 < 10) and (factr > 2.9 or factr < 1.1):
                                        type = 4
                                    elif ecc > 0.999 and factr < 3.14 and bg.b_area > 58 and l1/25 > l2:
                                        type = 3
                                    elif ecc > 0.999 and cdrp < 0.92 and cdrp > 0.86 and mlength > 23:
                                        type = 2
                                    elif ecc > 0.992 and ecc < 0.999 and bg.b_area < 50 and abs(ratx) > 0.96 and abs(ratx) < 0.98 and abs(raty) > 0.96 and abs(ratx) < 0.98:
                                        type = 4
                                    elif cdrp > 0.75 and cdrd < 0.182 and ((bg.b_area > 28) or (bg.b_area < 28 and mlength > 17)):
                                        if (ecc > 0.9996 or r < 0.028) and cdrp < 0.9 and mlength < 30 and bg.b_area < 62:#Worm catchers
                                            type = 2
                                        elif ecc > 0.99 and l1 > 400 and l1 < 600 and l2 > 60 and factr > 3.4:
                                            type = 2
                                        elif ecc < 0.99 and ecc > 0.975 and mlength < 17 and l1 < 16 and l2 > 2 and r > 0.2:
                                            type = 2
                                        elif ecc > 0.993 and factr < 3. and mlength < 40 and mlength > 28 and cdrp > 0.9 and cdrp < 0.94 and bg.b_area < 50:
                                            type = 2
                                        elif ecc > 0.993 and factr < 4 and factr > 3.5 and mlength < 25 and mlength > 17 and r < 0.12:
                                            type = 2
                                        elif ((factr < 3.76 and ecc > 0.99 and cdrd < 0.06 and r < 0.13 and (bg.b_area > 60. or mlength > 10.) and max(abs(ratx), abs(raty)) > 0.935) and (abs(ratx) > 0.9 and abs(raty) > 0.86) or
                                        factr < 4.1 and bg.b_area > 30 and cdrd < 0.059 and mlength < 16):
                                            type = 3
                                        elif (factr < 4.16 and cdrp > 0.74 and cdrd < 0.012 and bg.b_area < 50 and mlength < 20 and l1 < 23. and l1 > 12. and l2 < 3):
                                            type = 3
                                        else:
                                            type = 2
                                    elif cdrp > 0.75 and cdrd < 0.05 and bg.b_area > 30:
                                        type = 6
                                    elif cdrp < 0.6 and cdrp > 0.45 and cdrd < 0.5 and cdrd > 0.2 and ecc > 0.92 and ecc < 0.999:
                                        type = 3
                                    else:
                                        type = 2
                print_type(iid, type)
                #append = get_chars(str(type), 1)
                #fnom = files[0][:-1] if files[0][-1] == '/' else files[0]
                #if (not iid[-2:] == append):
                #    os.rename(efile, fnom + '/' + iid + append + '.' + tail)
            else:
                print('As of now only image analysis at contour level 40 is supported, sorry.') and exit()
