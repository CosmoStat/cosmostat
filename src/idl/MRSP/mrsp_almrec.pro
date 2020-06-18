;+
; NAME:
;        mrsp_almrec
;
; PURPOSE:
;   Computes the inverse spherical harmonic transform of a POLARIZED TQU map using the HEALPix representation (NESTED data).
;
; CALLING:
;     mrsp_almrec, Trans, imag, pixel_window=pixel_window
;
; INTPUTS:
;     Trans -- IDL structures of the ALM transform (see mrsp_almtrans.pro).
;
; OUTPUTS:
;     Imag -- IDL 2D array of a POLARIZED TQU healpix map image reconstructed
;
; KEYWORDS:
;
;     pixel_window -- scalar: if set, the image is convolved by the healpix pixel window (only for Healpix map)
;
; EXTERNAL CALLS:
;       anafast (healpix software)
;       alm2map_cxx (healpix C++ software)
;
; EXAMPLE:
;       Compute the inverse spherical harmonix transform  
;               mrsp_almrec, AlmTrans, ImagRec 
;         
; HISTORY:
;	Written: Pierrick Abrial & Jean-Luc Starck, 2005
;	December, 2005 File creation
;--------------------------------------------------------------------------------------------------------
;


;===============================================================

;===============================================================

pro alm_cxxitrans_p, Alm, img, nside=nside, nlmax=nlmax, index=index, pixel_window=pixel_window, fast=fast

COMMON MR1ENV

if not keyword_set(fast) then begin
  fast = DEF_ALM_FAST
  Niter = DEF_ALM_NITER
end else fast = 1

if not keyword_set(pixel_window) then pixel_window = 0

ALMFitsFile = gettmpfilename()
FileName = gettmpfilename()
alm2fits,index,Alm,ALMFitsFile

    command3 = gettmpfilename()
    openw, com2,command3,/get_lun
    printf,com2,'infile='+ALMFitsFile
    printf,com2,'outfile=!'+FileName
;    printf,com2,'almsfile='+ALMFitsFile
    printf,com2,'nside='+string(nside)
;    printf,com2,'plmfile=!plm.fits'
    printf,com2,'nlmax='+string(nlmax)
    printf,com2,'fwhm_arcmin =0'
    printf,com2,'polarisation=true'
    if keyword_set(fast) then printf,com2,'double_precision=false' $
    else printf,com2,'double_precision=true'
    if not keyword_set(pixel_window) then printf,com2,'pixel_window=false' $
    else  printf,com2,'pixel_window=true'
    free_lun,com2
  
 OutFileStdOut=gettmpfilename()
 spawn,'alm2map_cxx '+command3 + ' > ' + OutFileStdOut
; print, nlmax
 ; read_fits_map,FileName, img
  img = mrsp_read(FileName, /noverb)
  delete, FileName
  delete, command3
  delete, OutFileStdOut
  delete, ALMFitsFile
end

;===============================================================

pro mrsp_almrec, alm, imag, pixel_window=pixel_window

COMMON MR1ENV

if N_PARAMS() LT 2  then begin 
        print, 'CALLING SEQUENCE: mrs_alm, alm, Imag, to_healpix = to_healpix, to_glesp=to_glesp,pin=pin,nx=nx,np=np'
        goto, DONE
        end


alm2 = alm

if alm2.tab eq 1 then begin 

	;if alm2.complex_alm eq 1 then begin
		tab2alm_pola, alm2.alm, alm_coef, complex=alm2.complex_alm
	;end else begin
	;	tab2alm_pola, alm2.alm, alm_coef
	;end
	
end else begin

	alm_coef = alm2.alm

end

;print,'pixeltype',alm2.pixeltype

  ;print,'inv transform alm-> healpix'
  npix = alm2.npix
  nside = alm2.nside
  lmax = alm2.lmax
  
  if alm2.norm EQ 1 then alm_coef =  alm_coef / alm2.NormVal;	float(nside)
  
  alm_cxxitrans_p, alm_coef, imag, nside=nside, nlmax=alm2.lmax, index=alm2.index, pixel_window=pixel_window

DONE:

END



