#!/bin/bash
SCRIPT="runAnalysis.sh"
runNumber=1
while [ $runNumber -le 10 ]; do
    echo "Adding run number from file: $runNumber"

#make the script to submit
    (echo "#!/bin/bash"
echo "source /cvmfs/alice.cern.ch/etc/login.sh"
#VO_ALICE@ROOT::v6-24-06-3
echo "eval $(alienv printenv VO_ALICE@AliPhysics::vAN-20210910-1)"
echo "which aliroot || exit 1"
echo "mkdir -p /user/jlomker/project/AVFD/output/Set$runNumber"
#echo "cd /user/jlomker/project/AVFD/output/Set$runNumber"
#echo "mkdir -p ${TMPDIR}"
echo "cd ${TMPDIR}"
echo "pwd"
echo "if [ -f AnalysisResults.root ]"
echo "  then "
echo "rm -rf AnalysisResults.root"
echo "fi"
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
echo "exec aliroot -b -q runSingleCentrality.C'($runNumber)'"
    ) > $SCRIPT

qsub -q generic $SCRIPT 

let runNumber++

done
