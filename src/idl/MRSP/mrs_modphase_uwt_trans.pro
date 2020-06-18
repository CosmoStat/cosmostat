; NAME:
;        MRS_MODPHASE_UWT_TRANS
;
; PURPOSE:
;	Compute the UNdecimated wavelet transform on the sphere of a vector field (ex Q-U CMB data), 
;       using the healPix pixel representation (nested data representation). 
;       The discrete wavelet transform is applied successively on the  12 faces of the
;       Healpix image. On each face, the routine "mr_modphase_uwt_trans" is called. 
;       It computes first the modulus and the phase, and run an UNdecimated WT on each.
;       The output is a IDL structure.
;
; CALLING:
;     MRS_MODPHASE_UWT_TRANS, Imag, Trans, NbrScale=NbrScale
;
; INPUTS:
;     Imag -- IDL array of a healpix vector field  fltarr[*,2] : Input image be transformed 
;    
; OUTPUTS:
;     Trans -- IDL structures with the following fields:  
;         NBRSCALE  -- LONG: Number of scales of the wavelet transform
;         COEF      -- 3D IDL array [*,*,*,12,2] : Wavelet coefficients
;                        cube  containing all wavelet coefficients
;                        COEF[*,*, 0:NbrScale-1, 0:11, 0] = wavelet transform of the modulus.
;                        COEF[*,*, 0:NbrScale-1, 0:11, 1] = wavelet transform of the phase.
;	   Nx -- number of pixels on the side of the Healpix patch, nside
;	   Ny -- same as Nx			
;      Nside -- Nside value
;      Npix -- Number of pixels 
;
; KEYWORDS:
;         NBRSCALE  -- LONG: Number of scales of the wavelet transform
;
; EXTERNAL CALLS:
;         mr_modphase_uwt_trans
;
; EXAMPLE:
;       Compute the undecimated wavelet transform of a vector field I with five scales
;       The result is stored in WT
;               mrs_modphase_uwt_trans, Imag, WT, NbrScale=5
;               tvscl, WT.coef[*,*,j,f,0] ; plot the jth scale of the fth face wavelet transform (f = 0..11) of the modulus
;         
; HISTORY:
;	Written:  Jean-Luc Starck, May 2008
;-


pro mrs_modphase_uwt_trans, Imag, Trans,  NbrScale=NbrScale, modif=modif 

if N_PARAMS() LT 2  then begin 
        print, 'CALLING SEQUENCE: mrs_modphase_uwt_trans , Imag, Trans, NbrScale=NbrScale'
        goto, DONE
        end

if not keyword_set(NbrScale) then NbrScale = 4
if not keyword_set(modif) then modif = 0
 
npix = (size(imag))[1]
nside = npix2nside(npix)
Coef = fltarr(nside,nside, NbrScale, 12, 2)  


get_all_faces, Imag[*,0], CubeFace1
get_all_faces, Imag[*,1], CubeFace2

   
 ;loop over the 12 healpix base resolution pixels
 Face = fltarr(nside, nside, 2)
 for f=0,11 do begin
       Face[*,*,0] = CubeFace1[*,*,f]
       Face[*,*,1] = CubeFace2[*,*,f]
       mr_modphase_uwt_trans, Face, W, NbrScale=NbrScale
 	   Coef[*,*,*,f, 0] = W.MODCOEFF
       Coef[*,*,*,f, 1] = W.ANGCOEFF
 end
    
	   
Trans =  {nside:nside, npix:npix, NbrScale:NbrScale, Coef:Coef, Nx:nside, Ny:nside, modif:modif}
 

DONE:

END



