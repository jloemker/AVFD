##!/bin/bash
SCRIPT="runAnalysis.sh"
centID=2
while [ $centID -le 7 ]; do
echo "Opening the script for centrality ID: $centID"
​
#make the script to submit
    (echo "#!/bin/bash"
echo "source /cvmfs/alice.cern.ch/etc/login.sh"
​
#echo "eval $(alienv printenv AliPhysics/vAN-20210923_ROOT6-1)"
echo "eval $(alienv printenv VO_ALICE@AliPhysics::vAN-20210923_ROOT6-1)"
echo "which aliroot || exit 1"
echo "mkdir -p /data/alice/jlomker/AVFD/Centrality-$centID"
echo "cd /data/alice/jlomker/AVFD/Centrality-$centID"
echo "pwd"
echo "if [ ! -f convert_tree_splitFiles.C ]"
echo " then "
echo "ln -s /user/jlomker/project/AVFD/convert_tree_splitFiles.C ." 
echo "fi"
echo "if [ ! -f CalculateFlowCME.cxx ]"
echo " then "
echo "ln -s /user/jlomker/project/AVFD/CalculateFlowCME.cxx ."
echo "fi"
echo "if [ ! -f Particle.h ]"
echo " then "
echo "ln -s /user/jlomker/project/AVFD/Particle.h ."
echo "fi"
echo "if [ ! -f Event.h ]"
echo " then "
echo "ln -s /user/jlomker/project/AVFD/Event.h ."
echo "fi"
echo "if [ ! -f runSingleCentrality.C ]"
echo " then "
echo "ln -s /user/jlomker/project/AVFD/runSingleCentrality.C ."
echo "fi"
echo "exec aliroot -b -q convert_tree_splitFiles.C'($centID)' "
echo "exec aliroot -b -q CalulateFlowCME.cxx++ runSingleCentrality.C++'($centID)'"
echo "mv *.root /data/alice/user/jlomker/AVFD/Centrality-$centID/"
    )>$SCRIPT
​
qsub -q generic $SCRIPT 
​
let centID++
​
done
echo "success !"
