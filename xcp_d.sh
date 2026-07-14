#!/bin/bash
docker run -it \
-v /filepath/derivatives/fmriprep:/data:ro \
-v /filepath/xcp_d:/out \
-v /filepath/derivatives/codes/license.txt:/opt/freesurfer/license.txt \
-v /filepath/xcp_d/work:/work \
pennlinc/xcp_d:0.7.3 \
 /data /out \
 participant \
 --participant-label {subjID} \
 -p 36P \
 --lower-bpf 0.01 \
 --upper-bpf 0.08 \
 --smoothing 6 \
 -r 50 \
 -f 0.5 \
 --fs-license-file /opt/freesurfer/license.txt \
 --cifti \
 -w /work
