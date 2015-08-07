# DECOEventType_beta
ilhan's python event classification script
## Script info
### Usage

The script can be run in any python environment. The output is in the form [image id]:[event type].

To use, type into the terminal:
> /path/to/main.py path/to/img_folder --c 40

The script is not completely accurate and is being improved. Currently the average accuracy is 90 percent and the speed is approximately 1.864 seconds per jpeg image.

### Info

`main` and `main_txt` are folders containing the version of the script that runs with the folder as input and a textfile as input, respectively.

Each folder has the following contents (_ = corresponding name):

_.py : Original unedited script
_obf.py : Obfuscated and compressed script
_arc.pyz : Binary pyz file version of script

`hotspotsclass` contains the script to be run on the xandyCent.out file to identify hotspots.

All scripts output the predictions into the file "classifications.out" located in the folder you have `cd`ed into

### Accuracy

The script was run on a large sample of worms, tracks and spots:

Out of 497 worms it identified 471 correctly:<br>
94.5% accuracy<br>
<br>
Out of 123 tracks it identified 98 correctly (2 exempted for noise):<br>
79.7% accuracy<br>
<br>
Out of 88 spots it identified ~~76~~ 85 correctly (I looked through the wrong ones manually):<br>
96.6% accuracy<br>

Overall accuracy is then 90.27%

For every 539 events there are:<br>
363 worms,<br>
82 tracks and<br>
94 spots

List of [Supported image file](http://pillow.readthedocs.org/en/latest/handbook/image-file-formats.html) formats
