#!/bin/csh
#
source /cvmfs/cms.cern.ch/cmsset_default.csh
cd /uscmst1b_scratch/lpc1/old_scratch/lpceg/yurii/lovedeep/ggNtuples/Nov26ggNtuples/CMSSW_5_3_12/src/
cmsenv
cd ${_CONDOR_SCRATCH_DIR}

cp /uscmst1b_scratch/lpc1/old_scratch/lpceg/yurii/lovedeep/ggNtuples/Wgamma/skimmer_wg .
cp /uscmst1b_scratch/lpc1/old_scratch/lpceg/yurii/lovedeep/ggNtuples/Wgamma/Pileup*root .
cp /uscmst1b_scratch/lpc1/old_scratch/lpceg/yurii/lovedeep/ggNtuples/Wgamma/*.cc .
cp /uscmst1b_scratch/lpc1/old_scratch/lpceg/yurii/lovedeep/ggNtuples/Wgamma/*.h .
cp "/uscmst1b_scratch/lpc1/old_scratch/lpceg/yurii/lovedeep/ggNtuples/Wgamma/${2}" .

./skimmer_wg ${1} ${2}

echo "DONE!"




