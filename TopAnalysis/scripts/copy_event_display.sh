#! /bin/sh

rm -r ~/www/colour_flow/event_displays/*
cp ~/www/.htaccess ~/www/colour_flow/event_displays/
for dir in event_displays/*; 
do
    copy_dir=~/www/colour_flow/$dir
    mkdir $copy_dir
    cp $dir/*png $copy_dir
    cp ~/www/.htaccess $copy_dir
done