#! /bin/sh
echo "Begin file copy. file 1 $1 file 2 $2"

if [ ! -f $1 ];
then
    echo "NO FILE TO COPY"
    exit 1
fi
cp $1 $2
checksum1=`md5sum $1 | awk '{print $1;}'`
checksum2=`md5sum $2 | awk '{print $1;}'`
echo "checksum1 $checksum1 checksum2 $checksum2"
while [ "$checksum1" != "$checksum2" ];
do
    echo "Files differ"
    cp $1 $2
    echo "Copy repeated"
    checksum1=`md5sum $1 | awk '{print $1;}'`
    checksum2=`md5sum $2 | awk '{print $1;}'`
    echo "checksum1 $checksum1 checksum2 $checksum2"

done
echo "File copy done"
rm $1