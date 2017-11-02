#!/bin/bash
source movietime.functions
source movietime.settings
MODE=$1
DENPATH=$2

case $MODE in
	init)
		echo "Setting up the environment..."
		[ ! -f readwf-xy ] && Compile
		[ ! -f params ] && Compile
		[ ! -d Movie ] && makeDirStructure
		;;
	1)
		echo "Produce only densities..."
		for FILE in ${DENPATH}/density.*.dat
		do
			setFileID ${FILE}
			[ ! -f Movie/2D-densities/current.${ID}.dat ] && genCurrent ${FILE}
		done
		;;
	2)
		echo "Produce densities and images..."
		for FILE in ${DENPATH}/density.*.dat
		do
			setFileID ${FILE}
			[ ! -f Movie/2D-densities/current.${ID}.dat ] && genCurrent ${FILE}
			[ ! -f Movie/Params/param-${ID}.dat ] && genParams ${FILE}
			[ ! -f Movie/2D-images/denxz-${ID}.png ] && plotImage
			#[ ! -f Movie/1D-images/denz-${ID}.png ] && plot1DImage
		done
		;;
	3)
		echo "Produce densities, images and a movie..."
		for FILE in ${DENPATH}/density.*.dat
		do
			setFileID ${FILE}
			[ ! -f Movie/2D-densities/current.${ID}.dat ] && genCurrent ${FILE}
			[ ! -f Movie/Params/params.${ID}.py ] && genParams ${FILE}
			[ ! -f Movie/2D-images/den.${ID}.png ] && plotImagePy	
			#[ ! -f Movie/1D-images/denz-${ID}.png ] && plot1DImage
			echo
		done
		compileMovie
		;;
	4)
		echo "Produce a single density and image..."
		setFileID ${DENPATH}
		[ ! -f Movie/2D-densities/current.${ID}.dat ] && genCurrent ${DENPATH}
		if [[ $OVERWRITE == "true" ]]
		then
			genParams ${DENPATH}
			plotImagePy
			#plot1DImage
		else
			[ ! -f Movie/Params/params.${ID}.py ] && genParams ${DENPATH}
			[ ! -f Movie/2D-images/den.${ID}.png ] && plotImagePy
			#[ ! -f Movie/1D-images/denz-${ID}.png ] && plot1DImage
		fi
		;;
	5)
		echo "Produce only images..."
		for FILE in ${DENPATH}/density.*.dat
		do
			setFileID ${FILE}
			[ ! -f Movie/Params/param-${ID}.dat ] && genParams ${FILE}
			[ ! -f Movie/2D-images/denxz-${ID}.png ] && plotImage
			#[ ! -f Movie/1D-images/denz-${ID}.png ] && plot1DImage
		done
		;;		
	6)
		echo "Produce only a movie..."
		compileMovie
		;;
	cleanup)
		echo "Sanitising the environment..."
		rm -vf density params readwf-xy readwf-xz plot.pyc plotsettings.pyc djogger.dat
		cd src && make distclean
		;;
	*)
		clear
		echo
		echo "Please choose one of:"
		echo "====================="
		echo 
		echo "   ./movietime.sh  init  			: SETUP environment. This should always be ran FIRST"
		echo
		echo "   ./movietime.sh  1  /path/to/DENSITIES/	: produces 2D densities"
		echo "   ./movietime.sh  2  /path/to/DENSITIES/	: produces 2D densities + images"		
		echo "   ./movietime.sh  3  /path/to/DENSITIES/	: produces 2D densities + images + movie"		
		echo "   ./movietime.sh  4  /path/to/DENSITY.*.DAT	: produces only a single 2D density and image"
		echo "   ./movietime.sh  5  /path/to/DENSITIES/	: produces only the images (needs 2D densities)"
		echo "   ./movietime.sh  6  /path/to/IMAGES/		: produces only the movie (needs images)"		
		echo
		echo "   ./movietime.sh cleanup 			: SANITISE environment. Be CAREFUL, this will also"
		echo "						  delete the 'Movie' directory structure AND its CONTENTS."
		echo "						  You will have to run the SETUP-mode again afterwards."		
		echo
		;;
esac
