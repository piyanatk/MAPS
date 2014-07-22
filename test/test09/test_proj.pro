pro sph_trig,ra_ref,dec_ref,l,m
    sr = sin(ra_ref)
    sd = sin(dec_ref)
    cr = cos(ra_ref)
    cd = cos(dec_ref)
    
    theta = atan(l,m)
    dist  = sqrt(m^2 + l^2)
    dist  = asin(dist)
    x = cos(dist);
	y = sin(dist)*sin(theta);
	z = sin(dist)*cos(theta);

	; Now rotate x,y,z to desired ra,dec
		xp = -y*sr + cr*(x*cd - z*sd);
		yp = sr*(x*cd - z*sd) + y*cr;
		zp = x*sd + z*cd;

	; Fill appropriate grid point with (ra,dec)
	ra = atan(yp,xp);
	dec = atan(zp,sqrt(xp*xp+yp*yp));
    print,'ra: ',ra*180./!pi
    print,'dec: ',dec*180./!pi
end

pro interferometry,ra_ref,dec_ref,l,m
    sr = sin(ra_ref)
    sd = sin(dec_ref)
    cr = cos(ra_ref)
    cd = cos(dec_ref)

    dec=asin(m*cd + sqrt(1.-l^2-m^2)*sd)
    dra=atan(l/(sqrt(1.-l^2-m^2)*cd - m*sd))

    print,'ra: ',(dra+ra_ref)*180/!pi
    print,'dec: ',dec*180/!pi
end


