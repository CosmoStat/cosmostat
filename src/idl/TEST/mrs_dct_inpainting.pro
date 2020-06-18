function mrs_dct_inpainting, TabIn, Mask, opt=opt

f = h2f(TabIn)
fm = h2f(Mask)
fout = 0.*f

for ll=0,11 do begin
	
	print,'... Inpainting Face Nb ',  ll+1
	
	face = f[*,*,ll]
	vs = size(face)
	Nx = vs[1]
	facem = fm[*,*,ll]
	face = face * facem
  	
	ind = where(face ne 0.)
	;m = mean(face[ind])
	;face[ind] = face[ind] - m
	
	writefits,'face_in.fits',face;fits_write,'face_in.fits',face
 	
 	optMCA = ' -t5  -L2  -s0  -H  '
 	if keyword_set(opt) then optMCA = opt + optMCA 
 	
 	cmd = ' cb_mca  ' + optMCA  +  ' -B ' + strc(Nx) + ' face_in.fits face_out.fits'
	spawn, cmd
	
	tab = rim('face_out.fits')
 	fout[*,*,ll] =  tab 	 
endfor

InpMap = f2h(fout)
return, InpMap
end
