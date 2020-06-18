
 
;=========================================================================================

pro mrs_cmbsim, cl, NbrSim=NbrSim, files, noise=noise, nside=nside, map_no_noise=map_no_noise, Fwhm=Fwhm


if not keyword_set(nside) then nside = 256
if not keyword_set(NbrSim) then nb_real = 1 else nb_real = NbrSim
if not keyword_set(Fwhm) then Fwhm =  0.
npix = nside2npix (nside)
files = strarr(nb_real)
files_noise = strarr(nb_real)
cl_name= gettmpfilename() ; strcompress('cl_'+string(floor(1000000*randomu(seed)))+'.fits',/remove_all)
alm_name= gettmpfilename() ; strcompress('alm_'+string(floor(1000000*randomu(seed)))+'.fits',/remove_all)

lmax = (size(cl))[1]
if lmax ge 3L*nside then begin
   cl2 = cl[0:3L*nside]  
   lmax = 3L*nside
 end else cl2 = cl
cl2fits,cl2,cl_name
 
for i=0,nb_real-1 do begin
    files[i]=strcompress("real_"+string(i)+".fits",/remove_all)
    files_noise[i]=strcompress("real_noise_"+string(i)+".fits",/remove_all)
    command3 = gettmpfilename() ;strcompress('command3_'+string(floor(1000000*randomu(seed)))+'.dat',/remove_all)
    openw, com2,command3,/get_lun
    printf,com2,"infile="+cl_name
    printf,com2,'outfile=!'+alm_name
    printf,com2,'nlmax='+string(lmax)
    printf,com2,'fwhm_arcmin ='+string(Fwhm)
    printf,com2,'rand_seed='+string(floor(randomu(seed)*1000));-1'
    printf,com2,'polarisation=false'
    free_lun,com2
   spawn,'syn_alm_cxx '+command3;+ " >/dev/null"
   delete,command3

    command3 = gettmpfilename() ; strcompress('command3_'+string(floor(1000000*randomu(seed)))+'.dat',/remove_all)
    openw, com2,command3,/get_lun
    printf,com2,"infile="+alm_name
    printf,com2,'outfile=!'+files[i]
    printf,com2,'nlmax='+string(lmax)
    printf,com2,'nside='+string(nside)
    printf,com2,'polarisation=false'
    printf,com2,'pixel_window=false'
    free_lun,com2
   spawn,'alm2map_cxx '+command3;+ " >/dev/null"
   delete,command3
   delete,alm_name

   if keyword_set(noise) then begin 
                  read_fits_map,files[i],map
                  map = map + noise*randomn(seed,npix)
		  write_fits_map,files_noise[i],map,/ring
   endif
endfor

delete,cl_name  
end
