##!/bin/bash
SCRIPT="runAnalysis.sh"
dirID=0

#while [ $dirID -le 0 ]; do
echo "Directory ID: $dirID"
centID=1
while [ $centID -le 7 ]; do
pTiD=0
while [ $pTiD -le 5 ]; do
etaID=0
while [ $etaID -le 1 ]; do
echo "Opening the script for centrality ID: $centID , pT ID: $pTiD and eta ID: $etaID"
#make the script to submit
    (#echo "#!/bin/bash"
echo "source /cvmfs/alice.cern.ch/etc/login.sh"
echo "eval $(alienv printenv VO_ALICE@AliPhysics::vAN-20210923_ROOT6-1)"
#echo "eval $(alienv printenv VO_ALICE@ROOT::v6-24-06-1)"
#echo "which aliroot || exit 1"
echo "mkdir -p /data/alice/jlomker/AVFD/Centrality-$centID/dirID-$dirID"
echo "cd /data/alice/jlomker/AVFD/Centrality-$centID/dirID-$dirID"
echo "pwd"
echo "if [ -f *.root ]"
echo "  then "
echo "rm -rf *.root"
echo "fi"
echo "if [ ! -f convert_tree_splitFiles.C ]"
echo " then "
echo "ln -s /project/alice/users/jlomker/AVFD/convert_tree_splitFiles.C ." 
echo "fi"
echo "if [ ! -f CalculateFlowCME.h ]"
echo " then "
echo " ln -s /project/alice/users/jlomker/AVFD/CalculateFlowCME.h ."
echo "fi "
echo "if [ ! -f CalculateFlowCME.cxx ]"
echo " then "
echo "ln -s /project/alice/users/jlomker/AVFD/CalculateFlowCME.cxx ."
echo "fi"
echo "if [ ! -f Particle.h ]"
echo " then "
echo "ln -s /project/alice/users/jlomker/AVFD/Particle.h ."
echo "fi"
echo "if [ ! -f Event.h ]"
echo " then "
echo "ln -s /project/alice/users/jlomker/AVFD/Event.h ."
echo "fi"
echo "if [ ! -f runSingleCentrality.C ]"
echo " then "
echo "ln -s /project/alice/users/jlomker/AVFD/runSingleCentrality.C ."
echo "fi"
echo "exec root -b -q CalculateFlowCME.cxx++ runSingleCentrality.C++'($centID, $dirID, $pTiD, $etaID)'"
#echo "exec root -b -q convert_tree_splitFiles.C'($centID, $dirID)' CalculateFlowCME.cxx++ runSingleCentrality.C++'($centID,$dirID, $pTiD, $etaID)' "
    ) > $SCRIPT

qsub -q gpu $SCRIPT 

let etaID++

done 

let pTiD++

done

let centID++

done

#let dirID++

#done

echo "success !"
