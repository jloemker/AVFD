#!/bin/bash
SCRIPT="runAnalysis.sh"
#runNumber=1
#while [ $Cent_ID -le 8 ]; do
echo "Opening the script"

#make the script to submit
    (echo "#!/bin/bash"
echo "source /cvmfs/alice.cern.ch/etc/login.sh"

#echo "source /cmvfs/alice.cern.ch/bin/alienv"

#echo "eval $(alienv printenv VO_ALICE@AliPhysics::vAN-20210910-1)"
echo "eval $(alienv printenv AliPhysics/vAN-20210915_ROOT6-1)"
echo "which aliroot || exit 1"
#echo "mkdir -p /user/jlomker/project/AVFD/output"
#echo "cd /user/jlomker/project/AVFD/output"
#echo "mkdir -p ${TMPDIR}"
echo "cd ${TMPDIR}"
echo "pwd"
#echo "if [ -f AnalysisResults_$Cent_ID ]"
#echo "  then "
#echo "rm -rf AnalysisResults_$Cent_ID"
#echo "fi"
#echo "if [! -f runSingleCentrality.C]"
#echo "then "
#echo "ln -s /user/jlomker/project/AVFD/runSingleCentrality.C ."
#echo "fi"
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
echo "exec aliroot -b -q convert_tree_splitFiles.C"
echo "exec aliroot -b -q CalulateFlowCME.cxx"
echo "exec aliroot -b -q runSingleCentrality.C"
    ) > $SCRIPT
echo "success !"

qsub -q generic $SCRIPT 

#let runNumber++

#done
