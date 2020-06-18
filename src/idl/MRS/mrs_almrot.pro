;+
; NAME:
;        mrs_almrot
;
; PURPOSE:
;   Makes rotation of a map (or its alm decomposition) in the spherical harmonics domain, The outout is an array of Alm.  
;
; CALLING:
;     mrs_almrot, Map_or_Alm,  RotAlm, nsideRot=nsideRot, lmax=lmax,  alm=alm
;
; INPUTS:
;     Map_or_Alm -- IDL array of healpix map or IDL struture of an Alm decomposition : Input image to be analyzed 
;    
; OUTPUTS:
;     RotAlm -- IDL structures with the following fields: 
;					MapName -- string : map name 
;			        ALMRE -- IDL array(0:lmax, 0:lmax, 0:N-1) : real part of the rotated Alm,  N= nsideRot^2*12
;					ALMIM -- IDL array(0:lmax, 0:lmax, 0:N-1) : imaginary part of the rotated Alm,  N= nsideRot^2*12
;					THETAPHI -- IDL array(0:N-1, 5) 
;                                            THETAPHI[i,0] = Theta of of the ith rotation
;                                            THETAPHI[i,1] = Phi of of the ith rotation
;                                            THETAPHI[i,2] = x-coordinate of the vector corresponding to the ith rotation
;                                            THETAPHI[i,3] = y-coordinate of the vector corresponding to the ith rotation
;                                            THETAPHI[i,4] = z-coordinate of the vector corresponding to the ith rotation
;                                       the ith rotation corresponding to an angle (THETAPHI[i,0],  THETAPHI[i,1]) which corresponds alo to the angle
;                                       related to the ith pixel position  in the ring format.
;
; INPUT KEYWORDS:
;      lmax       :  integer -- lmax value in the analysis
;      nsideRot      :   resolution for angle rotations. Number of rotations is nsideRot^2*12
;
; EXTERNAL CALLS:
;       alm_t2i -- healpix routine
;       alm2fits -- healpix routine
;       mrs_alm_rotate -- C++ iSAP binary in $ISAP/bin
;
; EXAMPLE:
;       Compute the rotations correcponding to 256^2*12  angles, using lmax=5
;        mrs_almrot, Map,  RotAlm, nsideRot=256, lmax=5
;
;      Ditto, but using spherical harmonics for the input
;        mrs_almtrans, Map, Alm
;        mrs_almrot, Alm,  RotAlm, nsideRot=256, lmax=5, /Alm
;
;        mrs_almrot, Map_or_Alm,  RotAlm, nsideRot=nsideRot, lmax=lmax, norm=norm, alm=alm
;        mrs_almrot, Map_or_Alm,  RotAlm, nsideRot=nsideRot, lmax=lmax, norm=norm, alm=alm

 ;         
; HISTORY:
;       Anais Rassat & Jean-Luc Starck, 2013
;       May, 2013 File creation
;--------------------------------------------------------------------------------------------------------

pro mrs_almrot, Data,  RotAlm, nsideRot=nsideRot, lmax=lmax, norm=norm, alm=alm, new=new
COMMON C_PLANCK
COMMON MR1ENV

; alm_tab = alm.alm
; alm_t2i, alm_tab, index, almlist
; alm2fits, index, almlist, 'xx_test.fits'
;  spawn, "mrs_almrec -v xx_test.fits xx_rec.fits"
; r = mrs_read('xx_rec.fits')
RotAlm=0
map = Data
if keyword_set(alm) then  norm=0
if keyword_set(norm) then map = map / 1d6/2.725
if not  keyword_set(nsideRot) then begin
  if keyword_set (alm) then  nsideRot = Data.nside $
  else nsideRot = gnside(data)
endif

  FNin= gettmpfilename()
  FNout= gettmpfilename()
  cmd = BIN_ISAPCXX + '/mrs_alm_rotate '
   if keyword_set(new) then cmd = BIN_ISAPCXX + '/mrs_alm_rotate '
  if keyword_set(lmax) then cmd = cmd + ' -l ' + STRC(lmax)
  if keyword_set(nsideRot) then cmd = cmd + ' -n ' + STRC(nsideRot)
 if keyword_set(alm) then cmd = cmd + ' -a '
 
  cmd = cmd + ' ' +  FNin + ' ' + FNout
  if not keyword_set(alm) then mrs_write, FNin, map, /ring $
  else begin
     alm_tab = Data.alm
     alm_t2i, alm_tab, index, almlist
     alm2fits, index, almlist, FNin
  end
  
  ; print, cmd
  spawn, cmd
  FNre =  FNout + '_almre.fits'
  FNim =  FNout + '_almim.fits'
  FNAngle =  FNout + '_thetaphi.fits'
  
  Almre = readfits(FNre)
  Almim = readfits(FNim)
  ThetaPhi = readfits(FNAngle)
  vs = size(Almre)
  NRot = vs[3]
  ;help, vs
  ;print, nrot
  delete, FNin
  delete, FNre
  delete, FNim
  delete, FNAngle
 
RotAlm = {Info: 'RING indexing for rotations', NRot: NRot , NsideRot:NsideRot, Almre: Almre, Almim:Almim, ThetaPhi: ThetaPhi}
 
end
