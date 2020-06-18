pro mrsp_addnoise,map,mapdust,sigma
npix = (size(map))[1]
mapdust = map
for i= 0 ,2 do begin

mapdust(*,i) = map(*,i)+stddev(map(*,i))*sigma*randomn(seed,npix)

endfor

end


pro info_noise,dust,dust_noise,dust_rec
mrsp_addnoise,dust,dust_1,1
mrsp_addnoise,dust,dust_3,3
mrsp_addnoise,dust,dust_5,5

mrsp_wtfilter,dust_1,dust_teb1
mrsp_wtfilter,dust_3,dust_teb3
mrsp_wtfilter,dust_5,dust_teb5


mrsp_teb2tqu,dust_teb1,dust_tqu1
mrsp_teb2tqu,dust_teb3,dust_tqu3
mrsp_teb2tqu,dust_teb5,dust_tqu5


print,'bruit 1 sigma'
for i= 0 ,2 do begin
print,'echelle ',i
print,stddev(dust(*,i)),stddev(dust_1(*,i)-dust(*,i)),stddev(dust_tqu1(*,i)-dust(*,i))
endfor
print,'bruit 3 sigma'
for i= 0 ,2 do begin
print,'echelle ',i
print,stddev(dust(*,i)),stddev(dust_3(*,i)-dust(*,i)),stddev(dust_tqu3(*,i)-dust(*,i))
endfor
print,'bruit 5 sigma'
for i= 0 ,2 do begin
print,'echelle ',i
print,stddev(dust(*,i)),stddev(dust_5(*,i)-dust(*,i)),stddev(dust_tqu5(*,i)-dust(*,i))
endfor
save,filename='test_denoising_pola.sav',dust,dust_1,dust_3,dust_5,dust_teb1,dust_teb3,dust_teb5,dust_tqu1,$
dust_tqu3,dust_tqu5

print,'bruit 1 sigma'
for i= 0 ,2 do begin
print,'echelle ',i
mrs_wtfilter,dust_1(*,i),rec
print,stddev(dust(*,i)),stddev(dust_1(*,i)-dust(*,i)),stddev(rec-dust(*,i))
endfor
print,'bruit 3 sigma'
for i= 0 ,2 do begin
print,'echelle ',i
mrs_wtfilter,dust_3(*,i),rec
print,stddev(dust(*,i)),stddev(dust_3(*,i)-dust(*,i)),stddev(rec-dust(*,i))
endfor
print,'bruit 5 sigma'
for i= 0 ,2 do begin
print,'echelle ',i
mrs_wtfilter,dust_5(*,i),rec
print,stddev(dust(*,i)),stddev(dust_5(*,i)-dust(*,i)),stddev(rec-dust(*,i))
endfor


end
