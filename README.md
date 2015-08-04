# DECOEventType_beta
ilhan's python event classification script

## Script info

### Usage

The script can be run in any python environment. The output is in the form [image id]:[event type].

To use, type into the terminal:
> /path/to/main.py path/to/img_folder --c 40

The script is not completely accurate and is being improved. Currently the average accuracy is roughly 68 percent and the speed is about 2 seconds per jpeg image. Ambiguous events are not included in the accuracy calculation.

Details about accuracy:

The script was run on a large sample of worms and tracks...

Out of 32 worms it identified 30 correctly.<br>
Accuracy: 93.55%<br>
<br>
Out of 47 tracks it identified 18 correctly.<br>
Accuracy: 38.30%<br>
<br>
Out of 14 spots it identified 10 correctly.<br>
Accuracy: 71.43%<br>
<br>
The program identified all 3 ambiguous events correctly.<br>

There may have been some bias as the some of the events were rerun on the calibration sample. For these stats, Matt Meehan and I analyzed the events and compared our classification to the program's. The test was ["double blind"](https://explorable.com/double-blind-experiment) to ensure accurate stats.

For every 542 events there are (according to Heather's classifications):<br>
363 worms,<br>
82 tracks and<br>
94 spots

[Supported file formats](http://pillow.readthedocs.org/en/latest/handbook/image-file-formats.html): jpg jpeg bmp png eps gif im j2k j2p jpx msp pcx png ppm pgm pbm spi tiff webp xbm xv
