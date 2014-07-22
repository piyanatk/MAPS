#!/bin/bash
clear
echo "         Hi, $LOGNAME, welcome!           "
echo "                                          "
echo "                                          "
echo "##########################################"
echo "#--This is a script for MAPS simulation--#"
echo "#----      writen by Rusen Lu         ---#"
echo "##########################################"
echo "                                          "
echo "purpose:                                   "
echo "generating a uv-fits data given a sky image" 
echo "using the MAPS software                    "
echo "                                          "
#echo "usage: $0 YOUR_FITS_IMAGE"
#echo "the input image must have a suffix '.fits'  "
echo "                                           "

if test $# -lt 1 
then
echo "Two or more arguments are required"; exit
fi
help()
{
cat << HELP
USEAGE: $0 YOUR_FITS_IMAGE, the input image 
must has a suffix '.fits'
OPTIONS: -h help text
HELP
exit 0
}

[ -z "$1" ] && help
[ "$1" = "-h" ] && help

echo -e "Please enter the normalization factor for MAPS_im2uv\n"
echo -e "please make sure your normalizer is correct!\n"
echo -e "e.g., for 256*256 pixels, n=65536 \n"
read  normalizer
echo "you entered: $normalizer"

#normalizer=65536

echo "you are using normalizer = $normalizer"

sleep 1

MAPS_im2uv -i $1 -o  $(basename $1 .fits)_Visibility.dat -n $normalizer

sleep 2 
echo -e "now visgen is doing a noiseless simulation......\n"
sleep 2 
visgen -n $(basename $1 .fits) -s ALMA50 -A $SIM/array/sgra_array.txt -G  $(basename $1 .fits)_Visibility.dat  -V sgra_obs_spec4_256pixel -v -N
clear
sleep 1
echo -e "converting to uv-fits format.\n"
sleep 2
maps2uvfits $(basename $1 .fits).vis $(basename $1 .fits).uvfits -23.019278  -67.753167 5000  $SIM/array/sgra_array.txt
echo -e "Done!"
exit
