#!/bin/sh
data_dir="$CMSSW_BASE/src/Rivet/data_vbfz"
CMS_MGPythia="'Title={MG5+Py8 LO (CMS)}:LineColor=blue:ErrorBandColor=blue:ErrorBandOpacity=0.3:LineWidth=0.02'"
CMS_MGPythia_PATH="${data_dir}/vbfzjj_CMSPythia.yoda"
CMS_MGHerwigpp="'Title={MG5+Herwig++ LO (CMS)}:LineColor=green:ErrorBandColor=green:ErrorBandOpacity=0.3:LineWidth=0.02'"
CMS_Herwigpp_PATH="${data_dir}/vbfzjj_CMSHerwig.yoda"

comm_CMS_pythiaVherwig="rivet-mkhtml -n5 -o rivet_ZVBF_CMS_PythiaVHerwig -c style_ZVBF_LesHouchesStudy.plot $CMS_MGPythia_PATH:$CMS_MGPythia $CMS_Herwigpp_PATH:$CMS_MGHerwigpp"
echo $comm_CMS_pythiaVherwig
eval $comm_CMS_pythiaVherwig
