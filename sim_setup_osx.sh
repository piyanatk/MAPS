################################################################################
# sim_setup_osx.sh
#
# MAPS setup script template for Mac OS X
#
# 07/22/2014 : Piyanat Kittiwisit
#     Make this file based on sim_setup.sh
################################################################################

# I am not sure what these are for.
export SIM_CFLAGC="-O"
export SIM_CFLAGL=-O

################################################################################
# Set up environment variables needed to compile and run MAPS
################################################################################
if [ "$1" = "" ] ; then
    export MAPS_WCNAME='MAPS'
else
    export MAPS_WCNAME=$1
fi

# Change the next line to point to the base directory of the MAPS
export SIM=$HOME/src/$MAPS_WCNAME
# Uncomment the next line if you want terminal to tell you about MAPS setup
# echo MAPS directory is set to '$SIM='\'$SIM\'

# Other environment variables to make running MAPS easier
# The following environment variables seem to be used by internal codes.
export SIM_BIN=$SIM/bin
export SIM_LIB=$SIM/lib
export SIM_INC=$SIM/include
export VISDIR=$SIM/visdata

# The following environment variables are for convenient.
export TEXTDIR=$SIM/text
export ARRAYDIR=$SIM/array
export LAYOUTDIR=$SIM/stn_layout

################################################################################
# Set up environment variables for required libraries.
################################################################################
# Required libraries.
# Change the following lines to point to the locations of libraries.
# Mac users are recommended to use Homebrew to install these libraries.
export WCS_LIB=/usr/local/lib
export WCS_INC=/usr/local/include/wcslib
# Note: you might not need CBLAS to run MAPS.
export CBLAS_LIB=/usr/lib
export CBLAS_INC=/System/Library/Frameworks/Accelerate.framework/Versions/Current/Frameworks/vecLib.framework/Headers/
export CFITSIO_LIB=/usr/local/lib
export CFITSIO_INC=/usr/local/include
export NETCDF_LIB=/usr/local/lib
export NETCDF_INC=/usr/local/include

