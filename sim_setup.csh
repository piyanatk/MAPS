setenv SIM_CFLAGC "-O"
setenv SIM_CFLAGL -O


setenv HAVE_CBLAS FALSE 
setenv HAVE_WCSLIB TRUE


# you MUST change the line below to point to the base
# directory of the MAPS installation
setenv SIM $SOFT/MAPS

# you MUST change the lines below to point to the location of
# header files and libraries for wcslib and CBLAS

if ( $HAVE_WCSLIB == "TRUE" ) then
  echo "HAVE_WCSLIB is TRUE setting WCS_INC and WCS_LIB"
  setenv WCS_INC /usr/local/include/wcslib
  setenv WCS_LIB /usr/local/lib
else 
  echo "WARNING: HAVE_WCSLIB is set to FALSE not all GSM code will be compiled"
endif

if ( $HAVE_CBLAS == "TRUE" ) then
  echo "HAVE_CBLAS is TRUE setting  CBLAS_INC and CBLAS_LIB"
  setenv CBLAS_LIB 
  setenv CBLAS_INC 
else
  echo "WARNING: HAVE_CBLAS is set to FALSE - not all GSM code will be compiled"
endif

# you MIGHT have to modify these if cfitsio went into a place where
# your compiler will not automatically look
setenv CFITSIO_INC /usr/local/include
setenv CFITSIO_LIB /usr/local/lib

setenv SIM_LIB $SIM/lib
setenv SIM_INC $SIM/include
setenv VISDIR $SIM/visdata

setenv TEXTDIR $SIM/text
setenv ARRAYDIR $SIM/array	
setenv LAYOUTDIR $SIM/stn_layout

setenv GSMDIR $SIM/maps_makesky/gsm
setenv FILE408 $SIM/maps_makesky/scripts/408POLN.txt

# Avoid path bloat
if ("$?SIM_BIN" == "0") then
    setenv SIM_BIN ${SIM}/bin
    setenv PATH "${PATH}:${SIM_BIN}"
endif
