This script works in conjunction with main.py (a version of plotBlobs). 
It finds hotspots in images. A hotspot is a group of pixels that repetively are lit in a series of images.
They are in the same location in each image. Sometimes hotspots are the same across models of phones. 

Here are the instructions for this script.
1. Run the main.py script first! This prints x and y centers for each blob group to an output file that is creatively named 
xandyCent.out

2. Wait until main.py finishes

3. Run hotspotclass.py. It works if you are in the same folder as xandyCent.out

4. The output file, hotspotclass.out, should let you know if any of your images contain hotspots.
