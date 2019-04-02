#! /usr/bin/env bash

# Script for building the C components of SELPythonLibs in an optimal way
# Note that some libraries have outside dependencies that require 
# modification so CFLAGS0 and LDFLAGS0, depending on how you want to 
# satisfy them.


for dir in Stats_Lib CosmoCalcs Rotations BrokenPowerLaw Quadrature AssefTemplates MWReddening
do
	cd $dir

	echo ""
	echo "building in "$dir
	#setup environment variables
	if [[ `uname` == "Darwin" ]]; then
	    # Fink distribution flags
# 		export CFLAGS0="-I/sw/include"
# 		export LDFLAGS0="-L/sw/lib"
        # Mac Ports distribution flags
        export CFLAGS0="-I/opt/local/include"
        export LDFLAGS0="-L/opt/local/lib"
        
        IFS="."
	    versionarr=(`uname -r`)
	    unset IFS
	    if [[ ${versionarr[0]} -ge 18 ]]; then
	        #export CFLAGS0=$CFLAGS0" -I/Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX.sdk/usr/include"
	        export CFLAGS0="-isysroot /Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX.sdk "$CFLAGS0
	    fi
	    
	    export CFLAGS0=$CFLAGS0" -O3"
        
		if [[ `sysctl -a | grep sse` =~ ": 1" ]]; then
			export CFLAGS0=$CFLAGS0" -mfpmath=sse"
		fi

    	if [[ `sysctl -a | grep machdep.cpu.brand_string` =~ "Core(TM)2" ]]; then
			export CFLAGS0=$CFLAGS0" -mtune=core2 -march=core2"
    	elif [[ `sysctl -a | grep machdep.cpu.brand_string` =~ "Core(TM) i7" ]]; then
			export CFLAGS0=$CFLAGS0" -mtune=corei7 -march=corei7"
		else
			export CFLAGS0=$CFLAGS0" -mtune=native -march=native"
    	fi

    	if [[ `sysctl -a | grep "hw.optional.avx"` =~ ": 1" ]]; then
    		export CFLAGS0=$CFLAGS0" -mavx"
		fi

    	if [[ `sysctl -a | grep "hw.optional.fma"` =~ ": 1" ]]; then
			export CFLAGS0=$CFLAGS0" -mfma"
		fi

	elif [[ `uname` == "Linux" ]]; then
		export CFLAGS0="-std=c99 -O3 -fPIC -DUSE_RESTRICT"
		export LDFLAGS0="-L/usr/lib -lm"

		if [[ `cat /proc/cpuinfo | grep "model name"` =~ "Core(TM) i7" ]]; then
			export CFLAGS0=$CFLAGS0" -mtune=corei7 -march=corei7"
		else
			export CFLAGS0=$CFLAGS0" -mtune=native -march=native"
    	fi

		export CFLAGS0=$CFLAGS0" -ftree-loop-ivcanon"
		export CFLAGS0=$CFLAGS0" -funroll-loops"

		features=`cat /proc/cpuinfo | grep flags`
		if [[ $features =~ "sse" ]]; then
			export CFLAGS0=$CFLAGS0" -mfpmath=sse"
		fi

		if [[ $features =~ "avx" ]]; then
			export CFLAGS0=$CFLAGS0" -mavx"
		fi

    else
    	echo "Operating system not recognized"
	fi

	make clean
	make
	cd ../
done
