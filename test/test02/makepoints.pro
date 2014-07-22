
pro makepoints,
 dec = -26.5333 ; -26:32:00

 for i_az = 0,7 do begin
     for i_alt = 1,4 do begin
         az = i_az * 45.0
         alt = i_alt * 15.0
     endfor
 endfor

end
