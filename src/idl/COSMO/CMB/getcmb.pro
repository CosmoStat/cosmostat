; .r read_cambcl

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
   spawn,'echo $ISAP', path
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
; print, command
spawn,'syn_alm_cxx '+command + ' > ' + OutFileStdOut

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
lmax = long(nlmax)
x_sky = 0
y_sky = 0
np =0
nx = 0
TabNbrM = 0
norm1 = 0
normval =0d
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

 
alm = {PixelType:PixelType, tab: tab1, complex_alm: complex_alm, nside : long(nside), npix: long(nside2npix(nside)), ALM : Alm, norm: norm1, NormVal: double(NormVal),lmin:lmin,lmax:lmax, TabNbrM: TabNbrM, index:index , x_sky : x_sky , y_sky : y_sky , nx : nx , np : np}

end

;=========================================================================================


function getcl, wmap=wmap, ell=ell, nomonodip=nomonodip, verb=verb
if keyword_set(wmap) then Cl = readfits('$ISAP/data/def_cl.fits') $
else begin
      if keyword_set(Verb) then print, "PR1 Best fit planck_lens_lensedCls.dat"
      DIR = '$ISAP/param/cmb/'
      FN = 'planck_lens_lensedCls.dat'
      Cl = read_cambcl(dir+FN, ell=ell, nomonodip=nomonodip)
end
return, Cl
end

;=========================================================================================

; /monodip ==> remove the monopole and the dipole
; /LL1 ==> the Cl are multiplied bu l(l+1)/2PI, i.e. standard CAMB output
function read_pr1cl,  LL1=LL1, ell=ell, nomonodip=nomonodip
DIR = '$ISAP/param/cmb/'
FN = 'base_planck_lowl.bestfit_cl.dat'
Cl = read_cambcl(dir+FN, LL1=LL1, ell=ell, nomonodip=nomonodip)
return, Cl
end

;=========================================================================================

; /monodip ==> remove the monopole and the dipole
; /LL1 ==> the Cl are multiplied bu l(l+1)/2PI, i.e. standard CAMB output
function read_pr2cl,  LL1=LL1, ell=ell, nomonodip=nomonodip, pol=pol
DIR = '$ISAP/param/cmb/'
FN = 'COM_PowerSpect_CMB-base-plikHM-TT-lowTEB-minimum-theory_R2.02.txt'
Cl = read_cambclp (dir+FN, LL1=LL1, ell=ell, nomonodip=nomonodip, DATA_START=1)
if not keyword_set(pol) then Cl = Cl[*, 0]
return, Cl
end



;=========================================================================================

function getpolacmb, Cl=Cl, nside=nside, Fwhm=Fwhm, SigmaNoise=SigmaNoise, Cmb=Cmb, Alm=Alm, teb=teb, tab=tab, nomap=nomap
COMMON C_PLANCK

if not keyword_set(nside) then nside = 256
if not keyword_set(Fwhm) then Fwhm =  0.
if not keyword_set(Cl) then  Cl = read_pr2cl(/pol)

npix = nside2npix (nside)
cl_name= gettmpfilename() ;  strcompress('cl_'+string(floor(1000000*randomu(seed)))+'.fits',/remove_all)
lmax = (size(cl))[1]
if lmax GT P_LMAX then lmax = P_LMAX

if lmax ge 3L*nside then begin
   cl = cl[0:3L*nside, *] 
   lmax = 3L*nside
end

mrs_almsim, Alm, lmax, glesp=glesp, tab=tab, cl=cl, Fwhm=Fwhm, nside=nside, /pol
if keyword_set(nomap) then map = Alm $
else begin
     if not keyword_set(teb) then  mrsp_almrec, Alm, cmb $
     else begin
            alm_t = alm
            mrs_almrec, alm_t, cmb_t
            alm_e = alm
            alm_e.alm(*,*,0) = alm_e.alm(*,*,1)
            mrs_almrec,alm_e, e

            alm_b = alm
            alm_b.alm(*,*,0) = alm_b.alm(*,*,2)
            mrs_almrec,alm_b,b

            cmb = [[cmb_t], [e], [b]]
    end
    if keyword_set(SigmaNoise) then  map = cmb + SigmaNoise*randomn(seed,npix) else map = cmb
end

return, map
end


;===================================================

function getcmb, Cl=Cl, nside=nside, Fwhm=Fwhm, SigmaNoise=SigmaNoise, Cmb=Cmb, Alm=Alm, glesp=glesp, pol=pol, nomap=nomap, tab=tab, wmap=wmap, lmax=lmax, OLD=OLD, pr1=pr1
COMMON C_PLANCK

if keyword_set(pol) then map = getpolacmb(Cl=Cl, nside=nside, Fwhm=Fwhm, SigmaNoise=SigmaNoise, Cmb=Cmb, Alm=Alm, nomap=nomap, tab=tab)  $
else BEGIN

if not keyword_set(nside) then nside = 256
if not keyword_set(Fwhm) then Fwhm =  0.
if not keyword_set(Cl) then BEGIN
  ;  t = mrdfits('$ISAP/data/wmap_lcdm_bf_model_yr1_v1.fits', 1, Header) 
  ;  Cl = t.temperature
   if keyword_set(wmap) then Cl = readfits('$ISAP/data/def_cl.fits') $
   else begin
         if keyword_set(pr1) then begin
             ;DIR = '$ISAP/param/cmb/'
             ;FN = 'planck_lens_lensedCls.dat'
             ;Cl = read_cambcl(dir+FN)
             Cl = read_pr1cl()
         ; DIR = '$ISAP/param/cmb/PlanckParams/Planck/'
          ; FN = 'cl_lensed.dat'
         ; FN = 'taucl_lensed.dat'
          ; FN = 'taucl.dat'
           end else begin
              Cl = read_pr2cl()
          end
   end
   ; help, Cl
   END

npix = nside2npix (nside)
cl_name= gettmpfilename() ;  strcompress('cl_'+string(floor(1000000*randomu(seed)))+'.fits',/remove_all)
if not keyword_set(lmax) then lmax = (size(cl))[1]
if lmax GT P_LMAX then lmax = P_LMAX

if lmax ge 3L*nside then begin
   cl2 = cl[0:3L*nside] 
   lmax = 3L*nside
end else cl2 = cl
 
mrs_almsim, Alm, lmax, glesp=glesp, tab=tab, cl=cl2, Fwhm=Fwhm, nside=nside

if not keyword_set(nomap) then mrs_almrec, Alm, cmb $
else cmb = Alm 

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
END

return, map
end


;=========================================================================================

function getalmcmb, Cl=Cl, Fwhm=Fwhm,  pol=pol
z = getcmb(Cl=Cl, nside=nside, Fwhm=Fwhm,  pol=pol, /nomap)
return, z
end

;=========================================================================================
 
 
 
