;+
; NAME:
;        mrs_almtrans
;
; PURPOSE:
;   Computes the spherical harmonic transform, 
;   using the HEALPix representation (nested data
;   representation by default) or the GLESP representation. 
;
; CALLING:
;     mrs_almtrans, Imag, Trans, lmax=lmax, complex=complex, ring=ring, psp=psp, norm=norm, tab=tab
;
; INPUTS:
;     Imag -- IDL array of healpix map or GLESP structure: Input image to be transformed 
;    
; OUTPUTS:
;     Trans -- IDL structures with the following fields: 
;                      ALM: fltarray[*,2] = real and imaginary part of the ALM
;                           or  cfarr[*] = complex alm values if the keyword complex is set
;                           or  fltarray[NbrMaxM, NbrMaxL, 2] if the keyword tab is set
;                           or cfarr[NbrMaxM, NbrMaxL] of both keyword complex and tab are set.
;                      COMPLEX_ALM = 0 (default) ==> ALM contains the real and imaginary parts
;                                  = 1           ==> ALM contain an IDL complex array
;                                  = 2           ==> ALM contains the power spectrum and the phase
;                      PixelType: int = 0 for Healpix and 1 for GLESP
;                      tab: int = 0 for default ALM representation (i.e. 1D IDL array)
;                               = 1 for 2D representation (i.e. l for the first dimension and m for the second)
;                      nside : int = Healpix nside parameter, only used in Healpix representation
;                      lmax : int = maximum l value in the Spherical Harmonic Space.
;                      npix : int = Number of pixels of the input image (12*nside*nside for Healpix)
;                      complex_alm: int = 1 if the ALM values are in complex format
;                                         2 if the ALM values are in power spectrum and phase
;                      TabNbrM: IDL int array: set only of the /tab keyword is set.
;                      index: IDL int array: ALM pixel indices.
;
; KEYWORDS:
;      psp       : if set Trans.alm will contain the power spectrum and the phase instead of the real and imaginary parts 
;      complex   : if set Trans.alm will contain complex values instead of the real and imaginary parts 
;      Tab       : if set, ALM coefficients in Trans.alm are stored in a 2D array:
;                           Trans.alm[m,l]  where m = 0.. Trans.TabNbrM[l]-1  and l = 0..lmax-1
;      Lmax      : Number of spherical harmonics computed in the decomposition
;					(HEALPIX==> default is 3*nside, should be between 2*nside and 4*nside)
;					(GLESP==> default is: min([Imag.nx/2, Imag.np/4]
;      ring      : if set, the input data are in Heapix ring representation.
;      psp       : if set, the 
;      norm      : if set, a normalization is performed to the alm coefficient.
;
;
; EXTERNAL CALLS:
;       anafast (healpix software)
;       cl2map (glesp software)
;
; EXAMPLE:
;       Compute the spherical harmonix transform of an image. 
;        The result is stored in Output
;               mrs_trans, Imag, Output 
;         
; HISTORY:
;	Written: Pierrick Abrial & Jean-Luc Starck, 2005
;	December, 2005 File creation
;--------------------------------------------------------------------------------------------------------
 
;===============================================================

pro mrs_getalm_l2_l3, ima, l2, l3, l4, tit=tit, min=min, max=max, filename=filename
mrs_almtrans, ima, ai, lmax=10, /tab
al2 = ai
al2.alm[*]=0
al2.alm[2,*,*] = ai.alm[2,*,*]
al3 = ai
al3.alm[*]=0
al3.alm[3,*,*] = ai.alm[3,*,*]
al4 = ai
al4.alm[*]=0
al4.alm[4,*,*] = ai.alm[4,*,*]
mrs_almrec, al2, l2
mrs_almrec, al3, l3
mrs_almrec, al4, l4
T2=tit+ ': l=2'
T3=tit+ ': l=3'
T4=tit+ ': l=4'
FN2=filename+'_l2.png'
FN3=filename+'_l3.png'
FN4=filename+'_l4.png'
mrs_tv, l2, title=t2,  /healpix, min=min, max=max, png=FN2
mrs_tv, l3, title=t3,   /healpix, min=min, max=max, png=FN3
mrs_tv, l4, title=t4,  /healpix, min=min, max=max, png=FN4
end


;===============================================================

pro glesp_almt, data, tabAlm, NbrL=NbrL, index=index 

  if not keyword_set(NbrL)  then NbrL = min([(data.nx-1)/2,data.np/4])
  
  filename_tmp = gettmpfilename()  
  write_glesp, filename_tmp, data 
 
  command = 'cl2map  -s 1 -di+ -map '+  filename_tmp
  if keyword_set(NbrL) then command = command + " -lmax "+ string(NbrL)
  filename_alm = gettmpfilename() 
  command = command + " -ao " +  filename_alm
; print, command
  
  OutFileStdOut= gettmpfilename() 
  spawn, command + ' >& ' + OutFileStdOut  
  fits2alm, index, tabAlm, filename_alm
  delete, filename_tmp
  delete, filename_alm
  delete, 'cl2map.dat'
  delete, OutFileStdOut
end

;===============================================================

pro alm_trans, imag, Alm, nside=nside, nlmax=nlmax, index=index, ring=ring

npix = (size(imag))[1]
nside = npix2nside(npix)
if not keyword_set(nlmax) then begin
   nlmax = long( nside )  * 3l
   if nlmax GT 3000 then  nlmax = 3000L
end

FileName = gettmpfilename()  
ALMFitsFile = gettmpfilename() 
; FileName = 'tmp.fits'
; ALMFitsFile = 'alm.fits'				    
if not keyword_set(ring) then write_fits_map, FileName, imag, /nested $
else write_fits_map, FileName, imag, /ring

;analyse de l'image
; command = 'command.dat'
command = gettmpfilename() 
openw, com,command,/get_lun
printf,com,'infile='+FileName
printf,com,'outfile=!power.fits'
; printf,com,'plmfile=!plm.fits'
printf,com,'nlmax='+string(long(nlmax))
printf,com,'simul_type='+string(1)
;printf,com,'iter_order=5'
printf,com,'outfile_alms=!'+ALMFitsFile
free_lun,com

; print,"calcul des alm de l'image"
; print, command
OutFileStdOut=gettmpfilename() 
spawn,'anafast '+command + ' > ' + OutFileStdOut  

;recuperation des alm
fits2alm, index, Alm, ALMFitsFile
delete, FileName
delete, command
delete, OutFileStdOut
delete, ALMFitsFile
end

;===============================================================

pro alm_cxxtrans, imag, Alm, nside=nside, nlmax=nlmax, index=index, ring=ring, fast=fast
 COMMON C_PLANCK
COMMON MR1ENV

if not keyword_set(fast) then begin
  fast = DEF_ALM_FAST
  Niter = DEF_ALM_NITER
end else fast = 1

; print, 'CXX'
npix = (size(imag))[1]
nside = npix2nside(npix)
if not keyword_set(nlmax) then nlmax = long( nside )  * 3l

FileName = gettmpfilename() 
ALMFitsFile = gettmpfilename() 
if not keyword_set(ring) then write_fits_map, FileName, imag, /nested $
else write_fits_map, FileName, imag, /ring

;analyse de l'image
; command = 'command.dat'
spawn, "echo $HEALPIX/data", healp

command = gettmpfilename() 
openw, com,command,/get_lun
printf,com,'infile='+FileName
printf,com,'outfile=!power.fits'
; printf,com,'plmfile=!plm.fits'
printf,com,'nlmax='+string(long(nlmax))
printf,com,'simul_type='+string(1)
if not keyword_set(fast) then printf,com,'weighted=true' $  
else printf,com,'weighted=false'
printf,com,'healpix_data='+healp 
if not keyword_set(fast) then begin
  iter = 'iter_order=' + STRCOMPRESS(string(Niter), /REMOVE_ALL) 
 ;  printf,com, iter
  ; print, iter
end
if keyword_set(fast) then printf,com,'double_precision=false' $
else printf,com,'double_precision=true'
printf,com,'outfile_alms=!'+ALMFitsFile
printf,com,'polarisation=false'
free_lun,com

; print,"calcul des alm de l'image"
OutFileStdOut=gettmpfilename() 
; print, 'anafast_cxx '+command
; print, OutFileStdOut
spawn,'anafast_cxx '+command + ' > ' + OutFileStdOut
WAIT, P_WAIT

; print, nlmax

;recuperation des alm

fits2alm, index, Alm, ALMFitsFile
delete, FileName
delete, command
delete, OutFileStdOut
delete, ALMFitsFile
delete, 'power.fits'
end

;===============================================================
; l*l+l+m+1
function make_index, lmax
l = long(lmax)
m = l
Np = l*l+l+m+1l
i=0L
Tab = lonarr(Np)
for l=0l,lmax do begin
for m=0l,l do begin
Tab[i] =  l*l+l+m+1l
i=i+1
end
end
return, tab
end
;===============================================================

pro mrs_almtrans, Imag, out, lmax=lmax, ring=ring, tab=tab, complex=complex, psp=psp, norm=norm, fast=fast

COMMON MR1ENV
COMMON C_PLANCK

if N_PARAMS() LT 2  then begin 
        print, 'CALLING SEQUENCE: mrs_alm, Imag, out, lmax=lmax, ring=ring'
        goto, DONE
        end
out=0
PixelType=0
npix = 0L
nside = 0
complex_alm=0
TabNbrM = 0
x_sky = 0
y_sky = 0
nx  = 0
np  = 0
tab1=0
index=0l
if keyword_set(norm) then norm1 = 1 else norm1 = 0

if type_code(Imag) EQ 8 then begin
  PixelType=1
  x_sky = Imag.x_sky
  y_sky = Imag.y_sky 
  nx  = Imag.nx
  np = Imag.np
  vs = size(Imag.T_SKY)
  npix = vs[1]
  ;print,(Imag.nx-1)/2, (Imag.np)/4
  if not keyword_set(lmax)  then lmax = min([(Imag.nx-1)/2, (Imag.np)/4])
  glesp_almt, Imag, Alm, NbrL=lmax, index=index
 ;  info, alm[1:*,0]
  Alm = double(Alm)
  NormVal = double(npix) / lmax / 8.
  if norm1 EQ 1 then Alm = double(Alm) * NormVal
 ;  Alm[0,0]=mean(Imag.T_SKY)
end else begin
  npix = (size(imag))[1]
  nside = npix2nside(npix)
  if not keyword_set(lmax) then begin
     lmax = long( nside )  * 3l
     if lmax GT P_Lmax then  lmax = P_Lmax  
  end
  ; print, nside,  lmax

  if keyword_set(HealpixCXX) then begin
      if not keyword_set(ISAPCXX) then alm_cxxtrans, imag, Alm, nside=nside, nlmax=lmax, index=index, ring=ring, fast=fast  $
      else begin
        ;print, "NEW"
        FN1 = gettmpfilename()
        FN2 = gettmpfilename()
        mrs_write, FN1, Imag
        cmd = BIN_ISAPCXX + '/mrs_almtrans  -T ' 
	; cmd = BIN_ISAPCXX + '/mrs_almtrans2009  -T -p0'
        ;print, cmd
        if keyword_set(lmax) then cmd = cmd + ' -l ' + STRC(lmax)
        cmd = cmd + ' ' +  FN1 + ' ' + FN2
        ; print, cmd
        spawn, cmd
        A = readfits(FN2)
        tab2alm, A,  Alm
        index = make_index(lmax)
       ; Alm = mrdfits(FN2, 1)
       ; vs = size(alm.real)
       ; A = dblarr(vs[1], 2)
       ; A[*,0] = alm.real
       ; A[*,1] = alm.imag
       ; Alm = a
       ;  hs, alm
       delete, FN1
       delete, FN2
      end
  end else alm_trans, imag, Alm, nside=nside, nlmax=lmax, index=index, ring=ring
  ; help, Alm
  ; NormVal = double(nside)
  nel = float(N_ELEMENTS(Imag))
  NormVal =  sqrt(nel / (4.*!DPI))
  if norm1 EQ 1 then  Alm = double(Alm) *  NormVal
end
 
if keyword_set(complex) then begin
   complex_alm=1
   if not keyword_set(tab) then Alm = complex(Alm[*,0],Alm[*,1])
end else if keyword_set(psp) then begin
   complex_alm=2
   PSpectrum =  Alm[*,0]^2 + Alm[*,1]^2
   Phase = atan(Alm[*,1], Alm[*,0])
   Alm[*,0] = PSpectrum
   Alm[*,1] = Phase
end

; print, ' NbrL = ', lmax
; help, alm
if keyword_set(tab) then begin
    tab1=1
    alm2tab,Alm, TabALM, complex=complex, TabNbrM=TabNbrM, NbrL=lmax
    ALM =  TabALM
    end
out=0
lmin = 0
out = {PixelType:PixelType, tab: tab1, complex_alm: complex_alm, nside : nside, npix:npix, ALM : Alm, norm: norm1, NormVal: NormVal,lmin:lmin,lmax:lmax, TabNbrM: TabNbrM, index:index , x_sky : x_sky , y_sky : y_sky , nx : nx , np : np}

DONE:

END


