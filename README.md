# DECOEventType_beta
ilhan's python event classification script

## Script info

The script can be run in any python environemnt. The output is in the form [image id]:[event type].

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

For every 542 events there are (according to our classifications):<br>
363 worms,<br>
82 tracks and<br>
94 spots

As of now only .jpg and .jpeg images can be run using the script.
