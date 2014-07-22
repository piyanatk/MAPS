################################################################################
# sim_setup_u1204.sh
#
# MAPS setup script template for Ubuntu 12.04
#
# 07/19/2013 : Piyanat Kittiwisit
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

# Change the following lines to point to the base directory of the MAPS.
export SIM=$HOME/src/${MAPS_WCNAME}
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
# It seems that you do not need CBLAS to run MAPS.
export HAVE_CBLAS="FALSE"
# But WCSLIB is required.
export HAVE_WCSLIB="TRUE"
#export HAVE_WCSLIB="FALSE"

# Change the following lines to point to the location of header files and
# libraries for WCSLIB.
# The default location should work if you install WCSLIB through apt-get.
if [ "$HAVE_WCSLIB" = "TRUE" ] ; then
#  echo "HAVE_WCSLIB is set to TRUE"
  export WCS_INC=/usr/include/wcslib
  export WCS_LIB=/usr/lib
#else
#  echo "HAVE_WCSLIB is set to FALSE"
fi

# Change the following lines to point to the location of header files and
# libraries for CBLAS.
if [ "$HAVE_CBLAS" = "TRUE" ] ; then
#  echo "HAVE_CBLAS is set to TRUE"
  export CBLAS_LIB=
  export CBLAS_INC=
#else
#  echo "HAVE_CBLAS is set to FALSE"
fi

# Change the following lines to point to the location of header files and
# libraries for CFITSIO.
# The default location should work if you install CFITSIO through apt-get.
export CFITSIO_INC=/usr/include
export CFITSIO_LIB=/usr/lib/x86_64-linux-gnu
