###this script works in conjunction with main.py
###It analyzes the x and y coordinates of event candidates and then determines if the blob ("event") is a hotspot or an event.
###If the coordinates are repeated several times, then the blob is a hotspot. 
###It prints any images with hotspots to another file (in whatever folder you're in) called "hotspotsclass.out"


infile = open('xandyCent.out', 'r')
f = open('hotspotsclass.out', 'w')

#these are emtpy arrays for x and y centers
event_ID, xCent, yCent = [], [], []

#these are empty arrays that the script will add values into 
hotspot_xCent = []
xhotspot_img_index = []
hotspot_yCent = []
yhotspot_img_index = []
hotspots = []


#appending values into arrays
for line in infile:
    data = line.split(' ')
    event_ID.append(str(data[0]))
    xCent.append(float(data[1]))
    yCent.append(float(data[2]))

#gathering packages
import numpy as np
from collections import Counter
from collections import defaultdict

#round the coordinates. i.e. 732.987654321 would go down to 730. As hotspots tend to be in the same locatin (up to 5 pixel difference) this accounts for most hotspots.
xCent_rounded = np.around(xCent, decimals = -1)
yCent_rounded = np.around(yCent, decimals = -1)

#find duplicated xCent values w/ indexes: any repeated blob centers in each column (x) will be found
hotspot_xCent = defaultdict(list)
for i, item in enumerate(xCent_rounded):
	hotspot_xCent[item].append(i)
for k,v in hotspot_xCent.items():
	if len(v)>2:
		xhotspot_img_index.extend(v)
hotspot_xCent = sorted(xhotspot_img_index)

#find duplicated yCent values w/ indexes: any repated blob centers in each row will be found
hotspot_yCent = defaultdict(list)
for i, item in enumerate(yCent_rounded):
    hotspot_yCent[item].append(i)
for k,v in hotspot_yCent.items():
    if len(v)>2:
        yhotspot_img_index.extend(v)
hotspot_yCent = sorted(yhotspot_img_index)


#put them together: find intersections of repeated values. There will be a lot of blob centers at the same point if there is a hotspot.
hotspot_xyCent = hotspot_xCent + hotspot_yCent

#find duplicated coordinate values (hotspots) and their indexes.
h = Counter(hotspot_xyCent)

hotspots = [i for i in h if h[i] >= 3]

nonhotspots = [i for i in h if h[i] < 3]

#classification with event_ID

for i in range(0, len(hotspots)):
    print >>f, event_ID[hotspots[i]], "hotspot", (xCent[i], yCent[i])

print 'to see eventIDs with hotspots, look at hotspotclass.out file in your folder'
