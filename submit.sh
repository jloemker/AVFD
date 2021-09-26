##!/bin/bash
SCRIPT="runAnalysis.sh"
centID=0
while [ $centID -le 7 ]; do
echo "Opening the script for centrality ID: $centID"

#make the script to submit
    (echo "#!/bin/bash"
echo "source /cvmfs/alice.cern.ch/etc/login.sh"

echo "eval $(alienv printenv AliPhysics/vAN-20210915_ROOT6-1)"
echo "which aliroot || exit 1"
echo "cd /user/jlomker/project/AVFD"
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
echo "if [ ! -f Event.h]"
echo " then "
echo "ln -s /user/jlomker/project/AVFD/Event.h ."
echo "fi"
echo "if [ ! -f runSingleCentrality.C ]"
echo " then "
echo "ln -s /user/jlomker/project/AVFD/runSingleCentrality.C ."
echo "fi"
echo "exec aliroot .x convert_tree_splitFiles.C '($centID)'"
echo "exec aliroot .x CalulateFlowCME.cxx++"
echo "exec aliroot .x runSingleCentrality.C++ '($centID)'"
    ) > $SCRIPT

qsub -q generic $SCRIPT 

let centID++

done
echo "success !"
