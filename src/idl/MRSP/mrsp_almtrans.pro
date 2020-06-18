;+
; NAME:
;        mrsp_almtrans
;
; PURPOSE:
;   Computes the spherical harmonic transform of a POLARIZED TQU map using the HEALPix representation (NESTED data).
;
; CALLING:
;     mrsp_almtrans, Imag, Trans, lmax=lmax, tab=tab, complex=complex, norm=norm, fast=fast
;
; INPUTS:
;     Imag -- IDL array of healpix map or GLESP structure: Input image to be transformed 
;    
; OUTPUTS:
;     Trans -- IDL structures with the following fields: 
;                      ALM: fltarray[ *, 2, 3 ] = real and imaginary part of the ALM, ALM[*,*,0] is ALM T, ALM[*,*,1] is ALM E and ALM[*,*,2] is ALM B
;                           or cfarr[ *, 3 ] = complex alm values if the keyword complex is set
;                           or fltarray[ NbrMaxM, NbrMaxL, 2, 3 ] if the keyword tab is set
;                           or cfarr[ NbrMaxM, NbrMaxL, 3 ] of both keyword complex and tab are set.
;                      COMPLEX_ALM = 0 (default) ==> ALM contains the real and imaginary parts
;                                  = 1           ==> ALM contain an IDL complex array
;
;                      PixelType: int = 0 for Healpix (1 for GLESP but not used)
;                      tab: int = 0 for default ALM representation (i.e. 1D IDL array)
;                               = 1 for 2D representation (i.e. l for the first dimension and m for the second)
;                      nside : int = Healpix nside parameter
;                      lmax : int = maximum l value in the Spherical Harmonic Space.
;                      npix : int = Number of pixels of the input image (12*nside*nside for Healpix)
;
;                      TabNbrM: IDL int array: set only of the /tab keyword is set.
;                      index: IDL int array: ALM pixel indices.
;					   NormVal: float, normalization value aplied to the alm coefficients (if keyword /norm used)
;					   norm: int = 0 if no normalization has been aplied, else = 1
;
; KEYWORDS:
;      complex   : if set Trans.alm will contain complex values instead of the real and imaginary parts 
;      Tab       : if set, ALM coefficients in Trans.alm are stored in a 2D array:
;                           Trans.alm[m,l]  where m = 0.. Trans.TabNbrM[l]-1  and l = 0..lmax-1
;      Lmax      : Number of spherical harmonics computed in the decomposition HEALPIX==> default is 3*nside, should be between 2*nside and 4*nside
;      norm      : if set, a normalization is performed to the alm coefficient.
;
;
; EXTERNAL CALLS:
;       anafast (healpix software)
;
;
; EXAMPLE:
;       Compute the spherical harmonics transform of an image. 
;        The result is stored in Output
;               mrsp_almtrans, Imag, Output 
;         
; HISTORY:
;	Written: Pierrick Abrial & Jean-Luc Starck, 2005
;	December, 2005 File creation
;--------------------------------------------------------------------------------------------------------

;===============================================================

pro alm_cxxtrans_p, imag, Alm, nside=nside, nlmax=nlmax, index=index, fast=fast

COMMON MR1ENV

if not keyword_set(fast) then begin
  fast = DEF_ALM_FAST
  Niter = DEF_ALM_NITER
end else fast = 1

npix = (size(imag))[1]
nside = npix2nside(npix)
if not keyword_set(nlmax) then nlmax = long( nside )  * 3l

FileName = gettmpfilename() 
ALMFitsFile = gettmpfilename() 
mrsp_write,filename,imag

;analyse de l'image
spawn, "echo $HEALPIX/data", healp

command = gettmpfilename() 
openw, com,command,/get_lun
printf,com,'infile='+FileName
printf,com,'outfile=!power.fits'
; printf,com,'plmfile=!plm.fits'
printf,com,'nlmax='+string(nlmax)
printf,com,'simul_type='+string(1)
if not keyword_set(fast) then printf,com,'weighted=true' $  
else printf,com,'weighted=false'
printf,com,'polarisation=true'
printf,com,'healpix_data='+healp 
if not keyword_set(fast) then begin
  iter = 'iter_order=' + STRCOMPRESS(string(Niter), /REMOVE_ALL) 
  printf,com, iter
  ; print, iter
end
if keyword_set(fast) then printf,com,'double_precision=false' $
else printf,com,'double_precision=true'
printf,com,'outfile_alms=!'+ALMFitsFile
free_lun,com

; print,"calcul des alm de l'image"
OutFileStdOut=gettmpfilename() 
spawn,'anafast_cxx '+command + ' > ' + OutFileStdOut

; print, nlmax

;recuperation des alm
fits2alm, index, Alm, ALMFitsFile,'ALL'
delete, FileName
delete, command
delete, OutFileStdOut
delete, ALMFitsFile
delete, 'power.fits'
end

;===============================================================

pro mrsp_almtrans, Imag, out, lmax=lmax, tab=tab, complex=complex, norm=norm, fast=fast

COMMON MR1ENV
COMMON C_PLANCK

if N_PARAMS() LT 2  then begin 
        print, 'CALLING SEQUENCE: mrsp_almtrans, Imag, out, lmax=lmax, tab=tab, complex=complex, norm=norm, fast=fast'
        goto, DONE
        end
out=0
PixelType=0
npix = 0L
nside = 0L
TabNbrM = 0
complex_alm = 0

tab1=0
if keyword_set(norm) then norm1 = 1 else norm1 = 0

  npix = (size(imag))[1]
  nside = npix2nside(npix)
  
  if not keyword_set(lmax) then begin
     lmax = long( nside )  * 3l
     if lmax GT P_Lmax then  lmax = P_Lmax  
  endif
  ; print, nside,  lmax

 alm_cxxtrans_p, imag, Alm, nside=nside, nlmax=lmax, index=index,  fast=fast
  
  ;else alm_trans, imag, Alm, nside=nside, nlmax=lmax, index=index, ring=ring
  ; NormVal = double(nside)
  
  lm = sqrt( double(nside)*double(nside)*12.)
  NormVal = double(sqrt( double(lm*(lm+1)/(4.*!DPI)))) ;I THINK THIS I A BUG, should be sqrt((lm^2)/(4d *!dpi))
  
  if norm1 EQ 1 then  Alm = double(Alm) * NormVal
 
; print, ' NbrL = ', lmax
; help, alm

if keyword_set(complex) then begin

   complex_alm = 1
   if not keyword_set(tab) then begin
   		
   		alm_t = complex(Alm[*,0,0],Alm[*,1,0])
   		alm_e = complex(Alm[*,0,1],Alm[*,1,1])
   		alm_b = complex(Alm[*,0,2],Alm[*,1,2])
   		Alm = [ [alm_t], [alm_e], [alm_b] ]
   		
   	end
   
end

if keyword_set(tab) then begin
    
    tab1=1
    alm_pola2tab,Alm, TabALM, complex=complex, TabNbrM=TabNbrM, NbrL=lmax
    ALM = TabALM
    
end

out=0
lmin = 0
x_sky = 0
y_sky = 0
nx  = 0
np  = 0
out = { PixelType:PixelType, tab: tab1, complex_alm: complex_alm, nside : nside, npix:npix, ALM : Alm, norm: norm1, NormVal: NormVal, lmin:lmin, lmax:lmax, TabNbrM: TabNbrM, index:index, x_sky : x_sky , y_sky : y_sky , nx : nx , np : np}

DONE:

END

