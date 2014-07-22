
The NETCDF files can be read into IDL fairly easily for inspection/testing/debugging.
$ idl

fh = ncdf_open('test_turb.netcdf',/nowrite)
fi = ncdf_inquire(fh)
print,fi
names = strarr(fi.ndims)
sizes = lonarr(fi.ndims)
ids   = intarr(fi.ndims)
NCDF_DIMINQ, fh, 0, tempn, temps
names[0] = tempn
sizes[0] = temps
NCDF_DIMINQ, fh, 1, tempn, temps
names[1] = tempn
sizes[1] = temps
NCDF_DIMINQ, fh, 2, tempn, temps
names[2] = tempn
sizes[2] = temps
NCDF_DIMINQ, fh, 3, tempn, temps
names[3] = tempn
sizes[3] = temps

NCDF_VARGET,fh,0,dat
