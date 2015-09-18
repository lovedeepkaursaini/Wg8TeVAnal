#!/bin/csh
#
/bin/rm -f condor_${2}
cat > condor_${2} << +EOF

universe = vanilla
Executable = /uscmst1b_scratch/lpc1/old_scratch/lpceg/yurii/lovedeep/ggNtuples/Wgamma/ANA/SubmitCondor/makeProduction.csh
Should_Transfer_Files = YES
When_To_Transfer_Output = ON_EXIT
Output = /uscmst1b_scratch/lpc1/old_scratch/lpceg/yurii/lovedeep/ggNtuples/Wgamma/ANA/SubmitCondor/${2}.out
Error = /uscmst1b_scratch/lpc1/old_scratch/lpceg/yurii/lovedeep/ggNtuples/Wgamma/ANA/SubmitCondor/${2}.err
Log = /uscmst1b_scratch/lpc1/old_scratch/lpceg/yurii/lovedeep/ggNtuples/Wgamma/ANA/SubmitCondor/${2}.log
Arguments = ${1} ${2}
Queue 1

+EOF

/opt/condor/bin/condor_submit condor_${2}
