#! /bin/sh

rm -r LSFJOB_*
rm core*
rm LSF_outputs/*
rm analysis_muminus/*
rm analysis_muplus/*
rm analysis_munoniso/*
rm analysis/*

rm analysis_muminus/Chunks/*
rm analysis_muplus/Chunks/*
rm analysis_munoniso/Chunks/*
rm analysis/Chunks/*

rm analysis_muminus/plots/*
rm analysis_muplus/plots/*
rm analysis_munoniso/plots/*
rm analysis/plots/*

rm analysis_muminus/FARM/*
rm analysis_muplus/FARM/*
rm analysis_munoniso/FARM/*
rm analysis/FARM/*

rm migration/*

#source /afs/cern.ch/project/eos/installation/user/etc/setup.sh

rm -r $EOS/analysis/*
rm -r $EOS/LSF/*

rm -r $EOS/ReReco2016/f016290/migration
rm -r $EOS/ReReco2016/f016290/analysis_muplus
mkdir $EOS/ReReco2016/f016290/migration
mkdir $EOS/ReReco2016/f016290/analysis_muplus