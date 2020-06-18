;+
; NAME:
;        mrsp_wttrans
;
; PURPOSE:
;	Computes the undecimated isotropic wavelet transform of POLARIZED maps on the sphere, 
;   using the HEALPix representation (NESTED data representation). The wavelet function 
;	is zonal and its spherical harmonics coefficients a_l0 follow a cubic box-spline profile.
;	If DifInSH is set, wavelet coefficients are derived in the Spherical Harmonic Space, 
;	otherwise (default) they are derived in the direct space.
;
;
; CALLING:
;
; 		mrsp_wttrans, Imag, out, NbrScale=NbrScale, lmax=lmax, DifInSH=DifInSH, MeyerWave=MeyerWave     
;
; INPUTS:
;     Imag -- IDL array of healpix map: Input image to be transformed 
;    
; OUTPUTS:
;     Trans -- IDL structures with the following fields:  
;                  NbrScale : int = number of scales 
;                     nside : int = Healpix nside parameter
;                      lmax : int = maximum l value in the Spherical Harmonic Space (Healpix)
;                      npix : int = Number of pixels of the input image (12*nside*nside)
;                      Coef : fltarr[ npix, NbrScale, 3 ] = wavelet transform of the data. Coef[*,*,0] = wavelet transform on T, Coef[*,*,1] = wavelet transform on E, Coef[*,*,2] = wavelet transform on B
;                             Coef[ *, 0, *] = wavelet coefficients of the finest scale (highest frequencies).
;                             Coef[ *, NbrScale-1, *] = coarsest scale (lowest frequencies). 
;
; KEYWORDS:
;      NbrScale  : Number of scales (default is 4)
;      Lmax      : Number of spherical harmonics computed in the decomposition
;					(default is 3*nside, should be between 2*nside and 4*nside)
;      DifInSH   : If set, compute the wavelet coefficients as the
;					difference between two resolution in the spherical harmonics representation.
;					Otherwise, the wavelet coefficients are computed as the difference between two resolutions
;					in the initial representation.
;	   MeyerWave : If set, use Meyer wavelets and set the keyword DifInSH
;
; EXTERNAL CALLS:
;       anafast (healpix software)
;   	synfast (healpix software)
;   	alm_product2 (idl)
;   	compute_g (idl)
;   	compute_h (idl)
;
; EXAMPLE:
;
;       Compute the multiresolution of an image I with default options
;        The result is stored in Output
;               mrsp_wttrans, Imag, Output, NbrScale=5
;         
; HISTORY:
;	Written: Pierrick Abrial & Jean-Luc Starck, 2004
;	December, 2004 File creation
;--------------------------------------------------------------------------------------------------------

pro mrsp_wttrans, Imag, out, NbrScale=NbrScale, lmax=lmax, DifInSH=DifInSH, MeyerWave=MeyerWave

if N_PARAMS() LT 2  then begin 
        print, 'CALLING SEQUENCE: mrsp_wttrans, Imag, out, NbrScale=NbrScale, lmax=lmax, DifInSH=DifInSH, MeyerWave=MeyerWave'
        goto, DONE
        end
	 
if not keyword_set(NbrScale) then NbrScale = 4
if NbrScale le 1 or NbrScale ge 20 then begin print,'Error: Number of scales should be between 2 and 20'
    	    	    	    	    goto, DONE
				    end


nside=0

npix = (size(imag))[1]
nside = npix2nside(npix)
if not keyword_set(lmax) then lmax = nside *3


if not keyword_set(DifInSH)  then DifInSH = 0
if not keyword_set(MeyerWave) then MeyerWave = 0 else DifInSH = 1


ech = 1.
Hscale = imag  ; -> pb !!!!!

mrsp_almtrans, imag, ALM_p, lmax=lmax 

ALM_HighResolImag = ALM_p.alm
index=ALM_p.index
ech=1.
nlmax2 =lmax
							
TabWavelet = fltarr( npix, NbrScale, 3 )

for comp=0, 2 do begin 
	;boucle sur T, E , B 
	alm = alm_p
	hs, alm
	ech = 1.
	alm.alm(*,*,0) = alm_p.alm(*,*,Comp)
	alm.alm(*,*,1) *= 0
	alm.alm(*,*,2) *= 0
	ALM_HighResolImag = ALM.alm

	mrs_almrec,alm,Hscale ; correction ligne 85

	for j=0, NbrScale-2 do begin
		if keyword_set(MeyerWave) then begin
			hgmey, nlmax2, j, h, g, dif=dif 
		end else begin
			compute_g, lmax, ech, g
			compute_h, lmax, ech, h
		end
  
		alm_product2, ALM_HighResolImag, h, alm_h 
   
		if keyword_set(DifInSH) then begin
			alm_product2, ALM_HighResolImag, g, alm_g  
			ALM.alm = alm_g
			mrs_almrec, ALM, WScale
			TabWavelet[*,j,comp] = WScale 
		end else begin
			ALM.alm = alm_h
			mrs_almrec, ALM, LScale
			TabWavelet[*,j,comp] = Hscale - LScale
			Hscale = LScale
		end
		
		ALM_HighResolImag = alm_h
		ech = ech*2
	endfor

	j=NbrScale-1
	if keyword_set(DifInSH) then begin
		ALM.alm = alm_h
		mrs_almrec, ALM, LScale
		TabWavelet[*,j,comp] = LScale 
	end else begin 
		TabWavelet[*,j,comp] = LScale 
	end
endfor ; fin de boucle sur T,E,B

pyrtrans = 0
if not keyword_set(MeyerWave) then MeyerWave = 0 else MeyerWave = 1



TabNorm=[0.85,0.12,0.046,0.0224606,0.011,0.006] 
out = { NbrScale : NbrScale, nside : nside, npix : npix, Coef : TabWavelet, lmax : lmax, MeyerWave : MeyerWave, $
       DifInSH : DifInSH, pyrtrans : pyrtrans, TabNorm : TabNorm}
DONE:

END

