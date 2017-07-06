#1 /bin/sh
read -p "Enter queue number: " queue

#sh remove_plots.sh
#EOS=/store/cmst3/user/psilva/LJets2015/5736a2c
#EOS=/store/cmst3/user/psilva/LJets2015/b18c191
#EOS=/store/cmst3/user/psilva/LJets2015/076fb7a
#EOS_STORE=/store/cmst3/user/psilva/LJets2015/8c1e7c9
EOS_STORE=/store/cmst3/group/top/ReReco2016/b312177

#python scripts/runLocalAnalysis.py -i $EOS -q $queue --runSysts -o analysis_muplus   --ch 13   --charge 1
python scripts/runLocalAnalysis.py -i $EOS_STORE -q $queue  -o /eos/user/v/vveckaln/analysis  --era era2016 -m TOPJetShape::RunTopJetShape --skipexisting --only MC13TeV_TTJets --systVar all --exactonly


#python scripts/runLocalAnalysis.py -i $EOS -q $queue --runSysts -o analysis_muminus  --ch 13   --charge -1
#python scripts/runLocalAnalysis.py -i $EOS -q $queue            -o analysis_munoniso --ch 1300
