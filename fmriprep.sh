#!/bin/bash
docker run --security-opt seccomp=unconfined --rm -it \
-v /filepath:/data:ro \
-v /filepath/derivatives/fmriprep:/out \
-v /filepath/derivatives/codes/license.txt:/opt/freesurfer/license.txt \
-v /filepath/derivatives/fmriprep/work:/work \
nipreps/fmriprep:23.2.1 \
 /data /out \
 participant \
 --participant-label {subjID} \
 --use-aroma \
 --use-syn-sdc \
 --cifti-output \
 --dummy-scans 0 \
 --fs-license-file /opt/freesurfer/license.txt \
 --output-spaces MNI152NLin6Asym:res-2 T1w \
 -w /work
