# DECOEventType_beta
ilhan's python event classification script

## Script info

The script can be run in any python environment. The output is in the form [image id]:[event type].

The script is not completely accurate and is being improved. Currently the average accuracy is roughly 68 percent and the speed is about 2 seconds per jpeg image. Ambiguous events are not included in the accuracy calculation.

Details about accuracy:

The script was run on a large sample of worms and tracks (not spots intentionally)...

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

There may have been some bias as the some of the events were rerun on the calibration sample (which is becoming larger and larger...), and so the script will be rerun on a randomized sample. I am fairly certain that with some recalibration and retesting the script will have the same or higher accuracy on the random images as well. For these stats, Matt Meehan and I analyzed the events and compared our classification to the program's. The test was ["double blind"](https://explorable.com/double-blind-experiment) (Meaning the output was hidden until we agreed on the type) to ensure accurate stats.

For every 542 events there are (according to Heather's classifications):<br>
363 worms,<br>
82 tracks and<br>
94 spots

For now only .jpg and .jpeg images can be run using the script.
