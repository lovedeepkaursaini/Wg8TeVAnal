#!/bin/csh
#
source /cvmfs/cms.cern.ch/cmsset_default.csh
cd /uscmst1b_scratch/lpc1/old_scratch/lpceg/yurii/lovedeep/ggNtuples/Nov26ggNtuples/CMSSW_5_3_12/src/
cmsenv
cd ${_CONDOR_SCRATCH_DIR}

cp /uscmst1b_scratch/lpc1/old_scratch/lpceg/yurii/lovedeep/ggNtuples/Wgamma/ANA/run .
cp /uscmst1b_scratch/lpc1/old_scratch/lpceg/yurii/lovedeep/ggNtuples/Wgamma/ANA/*.cc .
cp /uscmst1b_scratch/lpc1/old_scratch/lpceg/yurii/lovedeep/ggNtuples/Wgamma/ANA/analyze.h .

./run ${1}

echo "DONE!"




