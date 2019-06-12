#!/bin/sh
data_dir="$CMSSW_BASE/src/Rivet/vbfz_data"
CMS_MGPythia="'Title={MG5+Py8 LO (CMS)}:LineColor=blue:ErrorBandColor=blue:ErrorBands=1:ErrorBandOpacity=0.3:LineWidth=0.02'"
CMS_MGPythia_PATH="${data_dir}/vbfzjj_pythia.yoda"
CMS_MGHerwigpp="'Title={MG5+Herwig++ LO (CMS)}:LineColor=green:ErrorBandColor=green:ErrorBands=1:ErrorBandOpacity=0.3:LineWidth=0.02'"
CMS_Herwigpp_PATH="${data_dir}/vbfzjj_herwigpp.yoda"

comm_CMS_pythiaVherwig="rivet-mkhtml -n5 -o rivet_VBFZ_CMS_PythiaVHerwig -c style_VBFZ_LesHouchesStudy.plot $CMS_MGPythia_PATH:$CMS_MGPythia $CMS_Herwigpp_PATH:$CMS_MGHerwigpp"
echo $comm_CMS_pythiaVherwig
eval $comm_CMS_pythiaVherwig
