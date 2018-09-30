#!/bin/bash
echo "source /data/highzgal/sourabh/bin/miniconda3/bin/activate astroconda; nohup jupyter-lab > /dev/null 2>&1 &" | ssh sourabh@ramon-4.spa.umn.edu:/bin/bash
# Obviously change the commands to start up jupyter lab etc.
ssh -N -f -L 8888:localhost:9999 sourabh@ramon-4.spa.umn.edu
#google-chrome "http://localhost:8888"
# "open" is for OSX but you can probably just invoke Chrome or something