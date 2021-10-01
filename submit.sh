##!/bin/bash
SCRIPT="runAnalysis.sh"
centID=2
while [ $centID -le 7 ]; do
echo "Opening the script for centrality ID: $centID"
#make the script to submit
    (echo "#!/bin/bash"
echo "source /cvmfs/alice.cern.ch/etc/login.sh"
#echo "eval $(alienv printenv AliPhysics/vAN-20210923_ROOT6-1)"
echo "eval $(alienv printenv VO_ALICE@AliPhysics::vAN-20210923_ROOT6-1)"
echo "which aliroot || exit 1"
echo "mkdir -p /data/alice/jlomker/AVFD/Centrality-$centID"
echo "cd /data/alice/jlomker/AVFD/Centrality-$centID"
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
#echo "exec aliroot -b -q convert_tree_splitFiles.C'($centID)' CalulateFlowCME.cxx++ runSingleCentrality.C++'($centID)'"
echo "exec aliroot -b -q convert_tree_splitFiles.C'($centID)' " 
    ) > $SCRIPT

qsub -q gpu $SCRIPT 

let centID++

done

echo "success !"
