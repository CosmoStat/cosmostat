function mrs_getpolacmb, Cl=Cl, nside=nside, Fwhm=Fwhm, SigmaNoise=SigmaNoise, Cmb=Cmb, Alm=Alm

if not keyword_set(nside) then nside = 256
if not keyword_set(Fwhm) then Fwhm =  0.
if not keyword_set(Cl) then BEGIN
   t = mrdfits('$MRS/data/wmap_lcdm_bf_model_yr1_v1.fits', 1, Header) 
   
   TT = t.temperature
   EE = t.gradient
   BB = t.curl
   TE = t.g_t
   
   Cl = [ [TT], [EE], [BB], [TE] ]
   help, Cl
   END

   

npix = nside2npix (nside)
cl_name= gettmpfilename() ;  strcompress('cl_'+string(floor(1000000*randomu(seed)))+'.fits',/remove_all)
lmax = (size(cl))[1]

if lmax ge 3L*nside then begin
   cl2 = cl[0:3L*nside] 
   lmax = 3L*nside
end else cl2 = cl
 
mrs_almsim, Alm, lmax, glesp=glesp, tab=tab, cl=cl1, Fwhm=Fwhm, nside=nside, /pol
mrsp_almrec, Alm, cmb

if keyword_set(SigmaNoise) then  map = cmb + SigmaNoise*randomn(seed,npix) else map = cmb
  
return, map
end


;=========================================================================================

function mrs_getcmb, Cl=Cl, nside=nside, Fwhm=Fwhm, SigmaNoise=SigmaNoise, Cmb=Cmb, Alm=Alm, glesp=glesp

if not keyword_set(nside) then nside = 256
if not keyword_set(Fwhm) then Fwhm =  0.
if not keyword_set(Cl) then BEGIN
   t = mrdfits('$MRS/data/wmap_lcdm_bf_model_yr1_v1.fits', 1, Header) 
   Cl = t.temperature
   ; help, Cl
   END

   

npix = nside2npix (nside)
cl_name= gettmpfilename() ;  strcompress('cl_'+string(floor(1000000*randomu(seed)))+'.fits',/remove_all)
lmax = (size(cl))[1]

if lmax ge 3L*nside then begin
   cl2 = cl[0:3L*nside] 
   lmax = 3L*nside
end else cl2 = cl
 
mrs_almsim, Alm, lmax, glesp=glesp, tab=tab, cl=cl1, Fwhm=Fwhm, nside=nside
mrs_almrec, Alm, cmb
if not keyword_set(Glesp) then  begin 
  if keyword_set(SigmaNoise) then  map = cmb + SigmaNoise*randomn(seed,npix) $
  else map = cmb
end else begin
  if keyword_set(SigmaNoise) then  begin
    vs = size(cmb.t_sky)
    map = cmb
    map.t_sky = cmb.t_sky + SigmaNoise*randomn(seed,vs[1]) 
  end else map = cmb
  end
  
return, map
end

;=========================================================================================

pro mrs_almsim, alm, nlmax, glesp=glesp, tab=tab, cl=cl, Fwhm=Fwhm, nside=nside, fast=fast, pol=pol
COMMON MR1ENV

if not keyword_set(pol) then Spol='polarisation=false'  else Spol='polarisation=true' 

if not keyword_set(fast) then begin
  fast = DEF_ALM_FAST
  Niter = DEF_ALM_NITER
end else fast = 1

if not keyword_set(Fwhm) then Fwhm =  0.

ALMFitsFile = gettmpfilename()
if keyword_set(cl) then begin
   CLFitsfile =gettmpfilename()
   cl2fits,cl,clfitsfile
endif else begin
   ALMFitsFile = gettmpfilename()
   spawn,'echo $MRS',path
   clfitsfile = path(0)+'/data/wmap_lcdm_bf_model_yr1_v1.fits'
endelse
command = gettmpfilename() 

openw, com,command,/get_lun
printf,com,'infile='+clfitsFile
printf,com,'outfile=!'+ALMFitsFile
printf,com,'nlmax='+string(nlmax)
printf,com,'rand_seed='+string(long(randomu(seed)*10000))
printf,com,'fwhm_arcmin ='+string(Fwhm)
if keyword_set(fast) then printf,com,'double_precision=false' $
else printf,com,'double_precision=true'
printf,com,Spol

free_lun,com

; print,"calcul des alm de l'image"
OutFileStdOut=gettmpfilename() 
cmd = 'syn_alm_cxx '+command + ' > ' + OutFileStdOut
; print, cmd
spawn, cmd
; print, "OK"

if not keyword_set(glesp) then pixeltype = 0 else pixeltype = 1

if not keyword_set(pol) then begin
	fits2alm, index, Alm, ALMFitsFile
end else begin
	fits2alm, index, Alm, ALMFitsFile,'All'
end

if not keyword_set(nside) then nside = 2^ceil(alog(nlmax /2.)/alog(2.))
nx = nlmax *2.
np = nlmax *4.
tab1 = 0
complex_alm = 0
lmin = 0
lmax = nlmax
x_sky = 0
y_sky = 0
np =0
nx = 0
TabNbrM = 0
norm1 = 0
normval =0
delete,ALMFitsFile
delete,command

if keyword_set(tab) then begin
	if not keyword_set(pol) then begin
		alm2tab,Alm, TabALM, complex=complex, TabNbrM=TabNbrM, NbrL=lmax
	end else begin
		alm_pola2tab, Alm, TabAlm, complex=complex, TabNbrM=TabNbrM, NbrL=lmax
	end	
	;if keyword_set(tab) then begin
	alm  = tabalm
	tab1=1
endif

alm = {PixelType:PixelType, tab: tab1, complex_alm: complex_alm, nside : nside, npix:nside2npix(nside), ALM : Alm, norm: norm1, NormVal: NormVal,lmin:lmin,lmax:lmax, TabNbrM: TabNbrM, index:index , x_sky : x_sky , y_sky : y_sky , nx : nx , np : np}

end



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
