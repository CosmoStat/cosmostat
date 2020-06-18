;+
; NAME:
;        mrs_attrans
;
; PURPOSE:
;		Compute the isotropic wavelet transform on the sphere, using the healPix pixel 
;		representation (NESTED data representation) and using the "A TROU" algorithm. 
;       The wavelet transform is applied successively on the  12 faces of the Healpix image.
;       The output is a IDL structure.
;
; CALLING:
;     mrs_attrans, Imag, Trans, NbrScale=NbrScale, Opt=Opt, modif=modif, healpix=healpix
;
; INPUTS:
;     Imag -- IDL array of healpix map: Input image be transformed 
;    
; OUTPUTS:
;     ATTrans -- IDL structures with the following fields:  
;         NBRSCALE  -- LONG: Number of scales of the wavelet transform
;          COEF      -- 4D IDL array [*,*,*,12] : Wavelet coefficients
;                        cube  containing all wavelet coefficients
;                        COEF[x,y, j, f] = wavelet coeff at face f (f=0..11),
;                                          position x,y and scale j
;	   Nx -- number of pixels on the side of the Healpix patch, nside
;	   Ny -- same as Nx			
;
; KEYWORDS:
;			NBRSCALE  -- LONG: Number of scales of the wavelet transform, default is 4
;			Opt -- string: if package MR1 is installed, extra keyword used by mr_transform.pro
;			modif: if set, add extra smoothing with spline filtering
;			healpix: if set, change ATTrans.coef to a 2D array [ c, j ] by reordering wavelet coefficients at scale j as a Healpix NESTED map.
;
; EXTERNAL CALLS:
;
; EXAMPLE:
;       Compute the orthogonal wavelet transform of an image I with five scales
;       The result is stored in Output
;               mrs_attrans, Imag, WT, NbrScale=5
;               tvscl, WT.coef[*,*,0,f] ; plot the fth face of the first scale of the wavelet transform  
;         
; HISTORY:
;	Written:  Jean-Luc Starck, 2006
;	October, 2006 File creation
;-

;-----------------------------------------------------------------

pro mrs_attrans, Imag, ATTrans,  NbrScale=NbrScale, Opt=Opt, modif=modif, healpix=healpix

COMMON MR1ENV

if N_PARAMS() LT 2  then begin 
        print, 'CALLING SEQUENCE: mrs_attrans, Imag, ATTrans,  NbrScale=NbrScale, Opt=Opt, modif=modif, healpix=healpix'
        goto, DONE
        end

if not keyword_set(NbrScale) then NbrScale = 4
if not keyword_set(modif) then modif = 0
if not keyword_set(healpix) then healpix = 0 else healpix = 1

npix = (size(imag))[1]
nside = npix2nside(npix)
Coef = fltarr(nside,nside, NbrScale, 12)  

get_all_faces, Imag, CubeFace

TabNorm=[0.85, 0.12, 0.046, 0.0224606, 0.011, 0.006, 0.006, 0.006, 0.006, 0.006] 
if keyword_set(modif) then TabNorm=[0.943169, 0.230368, 0.100762, 0.0488354, 0.0240847, 0.0120780, 0.00802764, 0.00802764, 0.00802764, 0.00802764] $
else  TabNorm=[0.890994,0.200981,0.0858255,0.0414706,0.0204158,0.0101239, 0.00802764, 0.00802764, 0.00802764, 0.00802764]


if keyword_set(mr1ok) then begin				;if package mre is available
print, 	'using mre'
	bin=1
	MR_FileName = 'xx_face.mr'
	OptWT=' '
	if keyword_set(Opt) then OptWT = OptWT + ' ' + Opt
	OptWT =  ' -n ' + STRCOMPRESS(string(NbrScale), /REMOVE_ALL) + OptWT 
	
	;loop over the 12 healpix base resolution pixels
	for f=0,11 do begin 
		mr_transform, CubeFace[*,*,f], Trans, Opt=OptWT, MR_File_Name=MR_FileName
		if f eq 0 then h = readfits(MR_FileName, HeadTrans)
		Coef[*,*,*,f] = Trans 
	end

	NbrScale = FXPAR(HeadTrans, "NBR_PLAN")
	if keyword_set(healpix) then begin
	   Coef1 = fltarr(npix,  NbrScale)
	   for j=0, NbrScale-1 do begin
	     Scale = reform(Coef[*,*,j,*])
	     Coef1[*,j] = F2H(Scale)
	   end
	   ATTrans = {Healpix: Healpix, nside: nside, npix: npix, UseGlesp: 0, pyrtrans: 0, NbrScale : NbrScale, Coef : Coef1, Nx: nside, Ny: nside, HeadTrans:HeadTrans, bin : 1, TabNorm: TabNorm}
	end else  ATTrans = {Healpix:  Healpix, nside: nside, npix: npix,  UseGlesp: 0, pyrtrans: 0, NbrScale : NbrScale, Coef : Coef, Nx: nside, Ny: nside, HeadTrans:HeadTrans, bin : 1, TabNorm: TabNorm}
 	delete, MR_FileName

endif else begin								;otherwise use procedures in the mrs package

   bin=0
   
   ;loop over the 12 healpix base resolution pixels
   for f=0,11 do begin
 		; W = bwt01_lift(CubeFace[*,*,f], NbrScale-1)
		atwt2d, CubeFace[*,*,f], W, Nscale=NbrScale,  modif=modif, /nostat
		Coef[*, *, *, f] = W 
	end
   if keyword_set(healpix) then begin
	   Coef1 = fltarr(npix,  NbrScale)
	   for j=0, NbrScale-1 do begin
	     Scale = reform(Coef[*,*,j,*])
	     Coef1[*,j] = F2H(Scale)
	   end
	   ATTrans =  {Healpix: Healpix,  nside: nside, npix: npix, UseGlesp: 0, pyrtrans: 0, NbrScale : NbrScale, Coef : Coef1, Nx: nside, Ny: nside, bin:0, modif: modif, TabNorm: TabNorm}
	end else  ATTrans = {Healpix: Healpix,  nside: nside, npix: npix, UseGlesp: 0, pyrtrans: 0, NbrScale : NbrScale, Coef : Coef, Nx: nside, Ny: nside, bin:0, modif: modif, TabNorm: TabNorm}
endelse


DONE:

END





pro mrs_atrec, Trans, HealpixIma 

COMMON MR1ENV


if N_PARAMS() LT 2 then begin 
        print, 'CALL SEQUENCE: mrs_atrec, WT_Struct, result'
        goto, DONE
        end
        
CubeFace = fltarr(Trans.Nx, Trans.Ny, 12)
NameSig = 'xx_signal.fits'

if Trans.Healpix EQ 1 then BEGIN
  Coef1 = fltarr (Trans.nside,Trans.nside, Trans.NbrScale, 12)  
  for j=0, Trans.NbrScale-1 do begin
             Scale = H2F(Trans.Coef[*,j])
	     Coef1[*,*,j,*] = Scale
  end
END

for f=0,11 do begin						;begin loop over 12 healpix faces

  if Trans.bin EQ 1 then begin				;if package mre is available and was used in the direct transform
		
		if Trans.Healpix EQ 0 then Coef = Trans.Coef[*,*,*,f] $
		else Coef = Coef1[*,*,*,f]
		; mr_recons, Coef, Rec, Header=Trans.HeadTrans
		atrec2d,  Coef, Rec
		CubeFace[*,*,f] = Rec 
    
  endif else begin							;use mrs for reconstruction as it was used for the direct transform
  
		if Trans.Healpix EQ 0 then Coef = Trans.Coef[*,*,*,f] $
		else Coef = Coef1[*,*,*,f]
 		atrec2d,  Coef, Rec, Modif=Trans.Modif
  		CubeFace[*,*,f] = Rec 
  endelse
end										;end loop over 12 healpix faces


put_all_faces, CubeFace, HealpixIma

delete, NameSig
DONE:
end


