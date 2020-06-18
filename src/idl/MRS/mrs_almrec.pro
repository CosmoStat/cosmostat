;+
; NAME:
;        mrs_almrec
;
; PURPOSE:
;   Computes the inverse spherical harmonic transform, 
;   using the HEALPix representation (nested data
;   representation by default) or the GLESP representation. 
;
; CALLING:
;     mrs_almrec, Trans, imag, to_healpix=to_healpix, to_glesp=to_glesp, pin=pin, nx=nx, np=np, pixel_window=pixel_window
;
; INTPUTS:
;     Trans -- IDL structures of the ALM transform (see mrs_almrec.pro).
;
; OUTPUTS:
;     Imag -- IDL array of healpix map or GLESP structure: Input image to be transformed 
;
; KEYWORDS:
;      to_healpix -- scalar: if true, force the output to be in Healpix format
;      to_glesp   -- scalar: if true, force the output to be in GLESP format
;      nx -- int: number of rings in GLESP (should be equal to 4nside)
;                 nx must b larger than 2* lmax
;      np -- int: number of pixels for the central ring in GLESP (should be equal to 8nside)
;                 np must b larger than 4 * lmax
;      pin : -1 equal aera   (default)
;	      0 isa lat/lon  (Theta-Phi image)
;	      1 triangular...etc..
;
;     ! warning : nx must >= 2 * nside !
;                 np must >= 4 * nside !
;
;     pixel_window -- scalar: if set, the image is convolved by the healpix pixel window (only for Healpix map)
;
; EXTERNAL CALLS:
;       anafast (healpix software)
;       cl2map (glesp software)
;       alm2map_cxx (healpix C++ software)
;
; EXAMPLE:
;       Compute the inverse spherical harmonix transform  
;               mrs_almrec, AlmTrans, ImagRec 
;         
; HISTORY:
;	Written: Pierrick Abrial & Jean-Luc Starck, 2005
;	December, 2005 File creation
;--------------------------------------------------------------------------------------------------------
;
pro glesp_alm_itrans, TabAlm, Imag, complex=complex, NbrL=NbrL,lmin=lmin,index=index, nx = nx, np = np, grn=grn
; Used
  ; tab2alm, TabAlm, alm,complex=complex
   filename_alm = gettmpfilename()
   alm2fits, index, float(TabAlm), filename_alm
    ; manque le tableau index !
   ;if not keyword_set(NbrL) then NbrL = 255
   if keyword_set(NbrL) and keyword_set(nx) and (Nbrl  ge nx /2) then NbrL = nx / 2
   if not keyword_set(NbrL) and keyword_set(nx) then NbrL = nx / 2
   
   filename_tmp = gettmpfilename()  
  ; dirG = '/Users/starck/Main/sapdev_osx_24_april_2013_before_cmake/external/glesp/'
  ; command =  'cl2map -falm  '+filename_alm+'  -o ' + filename_tmp +  '  -di+  -lmax '+ STRC(NbrL)
   command =  'cl2map -fa  '+filename_alm+'  -o ' + filename_tmp  +  ' -di+  -lmax '+ STRC(NbrL)  
   if keyword_set(lmin) then command = command + ' -lmin '+ STRC(lmin)
   if keyword_set(nx) then command = command + ' -nx ' + STRC(nx)
   if keyword_set(np) then command = command + ' -np ' + STRC(np)
   if not keyword_set(grn) then command = command + ' -gr a  '  $ ; quasi-equal Area pixelization 
   else command = command + ' -gr n  '   ; equal Number of pixels in each rings  

   OutFileStdOut=gettmpfilename()  
   ; print, command
   spawn, command+ ' >& ' + OutFileStdOut 
   read_glesp, filename_tmp,imag
   
   delete, filename_tmp
   delete, filename_alm
   delete, OutFileStdOut
   delete, 'alm.fits'
end



pro alm_itrans, Alm, img, nside=nside, nlmax=nlmax, index=index

; ALMFitsFile = 'alm.fits'				    
; FileName = 'tmp.fits'
ALMFitsFile = gettmpfilename()
FileName = gettmpfilename()
alm2fits,index,Alm,ALMFitsFile

    command3 = gettmpfilename()
    openw, com2,command3,/get_lun
    printf,com2,"infile=''"
    printf,com2,'outfile=!'+FileName
    printf,com2,'almsfile='+ALMFitsFile
    printf,com2,'nsmax='+string(nside)
;    printf,com2,'plmfile=!plm.fits'
    printf,com2,'nlmax='+string(nlmax)
    printf,com2,'fwhm_arcmin =0'
    printf,com2,'iseed=-1'
    free_lun,com2
  
 OutFileStdOut=gettmpfilename()
 spawn,'synfast '+command3 + ' > ' + OutFileStdOut
 ; read_fits_map,FileName, img
  img = mrs_read(FileName)
  delete, FileName
  delete, command3
  delete, OutFileStdOut
  delete, ALMFitsFile
end

;===============================================================

;===============================================================

pro alm_cxxitrans, Alm, img, nside=nside, nlmax=nlmax, index=index, pixel_window=pixel_window, fast=fast

COMMON C_PLANCK
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
    printf,com2,'polarisation=false'
    if keyword_set(fast) then printf,com2,'double_precision=false' $
    else printf,com2,'double_precision=true'
    if not keyword_set(pixel_window) then printf,com2,'pixel_window=false' $
    else  printf,com2,'pixel_window=true'
    free_lun,com2
  
 OutFileStdOut=gettmpfilename()
 spawn,'alm2map_cxx '+command3 + ' > ' + OutFileStdOut
 WAIT, P_WAIT

; print, nlmax
 ; read_fits_map,FileName, img
  img = mrs_read(FileName)
  delete, FileName
  delete, command3
  delete, OutFileStdOut
  delete, ALMFitsFile
end

;===============================================================

pro mrs_almrec, alm, imag,to_healpix = to_healpix, to_glesp=to_glesp, pin=pin, nx=nx, np=np, pixel_window=pixel_window 

COMMON MR1ENV

if N_PARAMS() LT 2  then begin 
        print, 'CALLING SEQUENCE: mrs_alm, alm, Imag, to_healpix = to_healpix, to_glesp=to_glesp,pin=pin,nx=nx,np=np'
        goto, DONE
        end


alm2 = alm

if keyword_set(to_healpix) and alm2.pixeltype eq 0 then begin 
   print, 'Error: to_healpix from healpix!!!'
   goto,done
endif

if keyword_set(to_glesp) and alm2.pixeltype eq 1 then begin
   print, 'Error: to_glesp from glesp!!'
   goto,done
endif


if alm2.tab eq 1 then begin 
   if alm2.complex_alm eq 1 then tab2alm, alm2.alm, alm_coef, complex=alm2.complex_alm $
   else tab2alm, alm2.alm, alm_coef
   if alm2.complex_alm eq 2 then begin
      Pspec = alm_coef[*,0]
      ind = where( Pspec LT 0, c)
      if c GT 0 then Pspec[ind]=0.
      alm_coef[*,0] = sqrt(Pspec) * cos (alm_coef[*,1])    
      alm_coef[*,1] = sqrt(Pspec) * sin (alm_coef[*,1])
    endif
endif else begin
    alm_coef = alm2.alm
    if alm2.complex_alm eq 1 then  alm_coef = [real_part(alm_coef),imaginary(alm_coef)]
    if alm2.complex_alm eq 2 then begin
      Pspec = alm_coef[*,0]
      ind = where( Pspec LT 0, c)
      if c GT 0 then Pspec[ind]=0.
      alm_coef[*,0] = sqrt(Pspec) * cos (alm_coef[*,1])    
      alm_coef[*,1] = sqrt(Pspec) * sin (alm_coef[*,1])
    endif
endelse

if keyword_set(to_healpix) then begin
  ; print,'to_healpix'
  alm2.pixeltype = 0
  if to_healpix GT 1 then  alm2.nside =to_healpix $
  else begin
    log = ceil(alog(alm2.lmax/3.)/alog(2.))
    alm2.nside = 2^log
  end
  alm2.npix = nside2npix(alm2.nside)
  endif

if keyword_set(to_glesp) then begin
  alm2.pixeltype = 1
  ;alm.nx = alm.nside*6
  ;alm.np = alm.nside*12
  if alm2.nx EQ 0 then alm2.nx = alm2.lmax * 2
  if alm2.np EQ 0 then alm2.np = alm2.lmax * 4
  endif

if keyword_set(nx) then alm2.nx= nx
if keyword_set(np) then alm2.np= np

if alm2.PIXELTYPE EQ 1 then begin
   MeanVal = alm_coef[0,0]
  if alm2.norm EQ 1 then alm_coef =  alm_coef  / float(alm2.npix) * alm2.lmax * 8.
  glesp_alm_itrans, alm_coef, Imag, complex=alm2.complex_alm, NbrL=alm2.lmax,lmin = alm2.lmin,index=alm2.index, nx = alm2.nx, np = alm2.np
 ;  Imag.T_sky[*] = Imag.T_sky[*] + MeanVal
endif else begin
   npix = alm2.npix
  nside = alm2.nside
  lmax = alm2.lmax
  if alm2.norm EQ 1 then alm_coef =  alm_coef / alm2.NormVal  ; float(nside)
  if keyword_set(HealpixCXX) then begin
    if not keyword_set(isapcxx) then alm_cxxitrans, alm_coef, imag, nside=nside, nlmax=alm2.lmax, index=alm2.index, pixel_window=pixel_window $
    else begin
        FN1 = gettmpfilename()
        FN2 = gettmpfilename()
        alm2tab,  alm_coef, A
        writefits, FN1, A
        cmd = BIN_ISAPCXX + '/mrs_almrec  -T -n  ' + strc(nside) + ' '  +  FN1 + ' ' + FN2
      ;  print, cmd
       spawn, cmd
        Imag = mrs_read(FN2)
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
   end else alm_itrans, alm_coef, imag, nside=nside, nlmax=alm2.lmax, index=alm2.index;, ring=ring
end


DONE:

END



