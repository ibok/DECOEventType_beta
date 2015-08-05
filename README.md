# DECOEventType_beta
ilhan's python event classification script<br><br>
## Script info

### Usage

The script can be run in any python environment. The output is in the form [image id]:[event type].

To use, type into the terminal:
> /path/to/main.py path/to/img_folder --c 40

The script is not completely accurate and is being improved. Currently the average accuracy is roughly 68 percent and the speed is about 2 seconds per jpeg image. Ambiguous events are not included in the accuracy calculation.

### Info

main.py is the main analysis script, to be run on a given folder.<br>
main_text.py runs on text files with the same syntax.<br>
main_obf.py is the same as main_text.py, but is obfuscated with [pyminifier](https://github.com/liftoff/pyminifier)<br>

All scripts output the predictions into the file "classifications.out" located in the folder you have `cd`ed into

### Accuracy

The script was run on a large sample of worms, tracks and spots:

Out of 518 worms it identified 471 correctly:<br>
90.9% accuracy<br>
<br>
Out of 129 tracks it identified 67 correctly:<br>
51.9% accuracy<br>
<br>
Out of 88 spots it identified 76 correctly:<br>
86.4% accuracy<br>

Overall accuracy is then 76.4%

For every 542 events there are:<br>
363 worms,<br>
82 tracks and<br>
94 spots

List of [Supported image file](http://pillow.readthedocs.org/en/latest/handbook/image-file-formats.html) formats
