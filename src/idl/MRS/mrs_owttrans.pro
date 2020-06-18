;+
; NAME:
;        mrs_owttrans
;
; PURPOSE:
;	Compute the (bi-) orthogonal wavelet transform on the sphere, 
;       using the healPix pixel representation (NESTED data representation). 
;       The wavelet transform is applied successively on the 12 faces of the
;       Healpix image.
;       The output is a IDL structure.
;
; CALLING:
;     mrs_owttrans, Imag, OWTTrans, NbrScale=NbrScale, Opt=Opt
;
; INPUTS:
;     Imag -- IDL array of healpix map: Input image be transformed 
;    
; OUTPUTS:
;     OWTTrans -- IDL structures with the following fields:  
;         NBRSCALE  -- LONG: Number of scales of the wavelet transform
;          COEF      -- 3D IDL array [*,*,12] : Wavelet coefficients
;                        cube  containing all wavelet coefficients
;                        COEF[*,*, f] = wavelet transform of face f (f=0..11).
;	   Nx -- number of pixels on the side of the Healpix patch, nside
;	   Ny -- same as Nx			
;
; KEYWORDS:
;         NBRSCALE  -- LONG: Number of scales of the wavelet transform
;		  Opt -- string: if package MR1 is installed, extra keyword used by mr_transform.pro
;
; EXTERNAL CALLS:
;         bwt01_lift (written by N. Aghanim and O. Forni)
;
; EXAMPLE:
;       Compute the orthogonal wavelet transform of an image I with five scales
;       The result is stored in Output
;               mrs_owttrans, Imag, WT, NbrScale=5
;               tvscl, WT.coef[*,*,f] ; plot the fth face wavelet transform (f = 0..11) 
;         
; HISTORY:
;	Written:  Jean-Luc Starck, 2005
;	February, 2005 File creation
;-
;-----------------------------------------------------------------

pro mrs_owttrans, Imag, OWTTrans, NbrScale=NbrScale, Opt=Opt

COMMON MR1ENV


if N_PARAMS() LT 2  then begin 
        print, 'CALLING SEQUENCE: Imag, OWTTrans,  NbrScale=NbrScale, Opt=Opt'
        goto, DONE
        end

npix = (size(imag))[1]
nside = npix2nside(npix)
Coef = fltarr(nside,nside, 12)
get_all_faces, Imag, CubeFace

		
if keyword_set(mr1ok) then begin				;if package mre is available
print, 	'using mre'
	bin=1
	MR_FileName = 'xx_face.mr'
	OptWT='-t14 -L'
	if keyword_set(Opt) then OptWT = OptWT + ' ' + Opt
	if keyword_set(NbrScale) then OptWT = OptWT + ' -n ' + STRCOMPRESS(string(NbrScale), /REMOVE_ALL) 
	
	;loop over the 12 healpix base resolution pixels
	for f=0,11 do begin 
		mr_transform, CubeFace[*,*,f], Trans, Opt=OptWT, MR_File_Name=MR_FileName
		if f eq 0 then h = readfits(MR_FileName, HeadTrans)
		Coef[*,*,f] = Trans 
	end

	NbrScale = FXPAR(HeadTrans, "NBR_PLAN")
	OWTTrans = {NbrScale : NbrScale, Coef : Coef, Nx: nside, Ny: nside, HeadTrans:HeadTrans, bin : 1}

	delete, MR_FileName

endif else begin								;otherwise use procedures in the mrs package

   bin=0
   if not keyword_set(NbrScale) then NbrScale = 4
   
   ;loop over the 12 healpix base resolution pixels
   for f=0,11 do begin
		;W = BWT01_DIRECT( CubeFace[*,*,f], NbrScale, MOMENT)
		W = bwt01_lift(CubeFace[*,*,f], NbrScale-1)
		Coef[*,*,f] = W 
	end
	
	OWTTrans = {NbrScale : NbrScale, Coef : Coef, Nx: nside, Ny: nside, bin:0}
endelse




DONE:

END
   
;====================================================================================================

