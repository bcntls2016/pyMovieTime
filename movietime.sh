#!/bin/bash
source movietime.functions
source movietime.settings
MODE=$1
NAME=$2
DENPATH=$3
OUTDIR=$4

case $MODE in
	init)
		echo "Setting up the environment..."
		[ ! -f readwf-xy ] && Compile
		[ ! -f params ] && Compile
		[ ! -d ${OUTDIR} ] && makeDirStructure
		;;
	1) ## LOOPED MODE
		echo "Produce only DENSITIES..."
		for FILE in ${DENPATH}/density.*.dat
		do
			setFileID ${FILE}
			[ ! -f ${OUTDIR}/2D-densities/current.${ID}.dat ] && genCurrent ${FILE}
		done
		;;
	2) ## LOOPED MODE
		echo "Produce DENSITIES and IMAGES..."
		for FILE in ${DENPATH}/density.*.dat
		do
			setFileID ${FILE}
			[ ! -f ${OUTDIR}/2D-densities/current.${ID}.dat ] && genCurrent ${FILE}
			[ ! -f ${OUTDIR}/Params/param-${ID}.dat ] && genParams ${FILE}
			[ ! -f ${OUTDIR}/2D-images/denxz-${ID}.png ] && plotImage
			#[ ! -f ${OUTDIR}/1D-images/denz-${ID}.png ] && plot1DImage
		done
		;;
	3) ## LOOPED MODE
		echo "Produce DENSITIES, IMAGES and a MOVIE..."
		for FILE in ${DENPATH}/density.*.dat
		do
			setFileID ${FILE}
			[ ! -f ${OUTDIR}/2D-densities/current.${ID}.dat ] && genCurrent ${FILE}
			[ ! -f ${OUTDIR}/Params/params.${ID}.py ] && genParams ${FILE}
			[ ! -f ${OUTDIR}/2D-images/den.${ID}.png ] && plotImagePy	
			#[ ! -f ${OUTDIR}/1D-images/denz-${ID}.png ] && plot1DImage
			echo
		done
		compileMovie
		;;
	4) ## LOOPED MODE
		echo "Produce IMAGES and a MOVIE..."
		for FILE in ${DENPATH}/density.*.dat
		do
			setFileID ${FILE}
			[ ! -f ${OUTDIR}/Params/param-${ID}.dat ] && genParams ${FILE}
			[ ! -f ${OUTDIR}/2D-images/denxz-${ID}.png ] && plotImage
			#[ ! -f ${OUTDIR}/1D-images/denz-${ID}.png ] && plot1DImage
		done
		compileMovie
		;;		
	5)
		echo "Produce a single density and image..."
		makeDirStructure
		setFileID ${DENPATH}/${NAME}
		[ ! -f ${OUTDIR}/2D-densities/current.${ID}.dat ] && genCurrent ${DENPATH}
		if [[ $OVERWRITE == "true" ]]
		then
			genParams ${DENPATH}
			plotImagePy
			#plot1DImage
		else
			[ ! -f ${OUTDIR}/Params/params.${ID}.py ] && genParams ${DENPATH}
			[ ! -f ${OUTDIR}/2D-images/den.${ID}.png ] && plotImagePy
			#[ ! -f ${OUTDIR}/1D-images/denz-${ID}.png ] && plot1DImage
		fi
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
