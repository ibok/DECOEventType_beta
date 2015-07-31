# DECOEventType_beta
python event classification script (unfinished)

## Script info

Run the script in the terminal or wherever else you prefer. The output is in the form [image id]:[event type].

The script is not completely accurate and is being improved. Currently the weighted accuracy is roughly 81.28% percent and the speed is about 2 seconds per jpeg image. Ambiguous events are not included in the accuracy calculation.

Details about accuracy:

The script was run on a large sample of worms and tracks...

Out of 32 worms it identified 30 correctly.<br>
Accuracy: 93.55%<br>
<br>
Out of 47 tracks it identified 18 correctly.<br>
Accuracy: 38.30%<br>
<br>
Out of 14 spots it identified 10 correctly.<br>
(A larger sample of pure spots might yield a higher accuracy)<br>
Accuracy: 71.43%<br>
<br>
The program identified all 3 ambiguous events correctly.<br>

For every 542 events there are, according to our classifications:<br>
363 worms,<br>
82 tracks and<br>
94 spots

The weighted score reflects this distribution. The formula for the weighted score is: (percent correct * frequency)[for each type] divided by / all frequencies added together.

Note: Only .jpg and .jpeg images (caps work too) can be run using the script.

[http://wipac.wisc.edu/deco](http://wipac.wisc.edu/deco)
