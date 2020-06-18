;+
; NAME:
;        mrsp_wtrec
;
; PURPOSE:
;	Compute the inverse wavelet transform of POLARIZED maps on the sphere.
;
;
; CALLING:
;
;     mrsp_wtrec, Trans, Rec, filter=filter        
;    
; INPUT:
;     Trans -- IDL structures with the following fields:  
;                  NbrScale : int = number of scales 
;                     nside : int = Healpix nside parameter
;                      npix : int = Number of pixels
;                      Coef : fltarr[ npix, NbrScale, 3 ] = wavelet transform of the data. Coef[*,*,0] = wavelet transform on T, Coef[*,*,1] = wavelet transform on E, Coef[*,*,2] = wavelet transform on B
;                             Coef[ *, 0, *] = wavelet coefficients of the finest scale (highest frequencies).
;                             Coef[ *, NbrScale-1, *] = coarsest scale (lowest frequencies).
;                       lmax: int = lmax parameter at the last scale
;
; OUTPUT:
;     Imag -- IDL array of healpix map: reconstructed image from the wavelet coefficients   
;
; KEYWORDS:
;      filter : Use filters for the reconstructions. If this keyword is not set, the reconstructed image
;               is obtained by a simple addition of all wavelet scales.
;
; EXTERNAL CALLS:
;       anafast (healpix software)
;   	synfast (healpix software)
;   	alm_product2 (idl)
;   	compute_g (idl)
;   	compute_h (idl)
;   	compute_gtilde (idl)
;   	compute_htilde (idl)
;
; EXAMPLE:
;       Compute the inverse wavelet transform:
;        The result is stored in Imag 
;               mrsp_wtrec, Trans, Imag 
;         
; HISTORY:
;	Written: Pierrick Abrial & Jean-Luc Starck, 2004
;	December, 2004 File creation
;------------------------------------------------------------------------------

pro mrsp_wtrec, Trans, Imag, filter=filter 

 
if N_PARAMS() LT 2 or N_PARAMS() GE 3 then begin 
        print, 'CALLING SEQUENCE: mrs_wtrec, Trans, Imag, filter=filter'
        goto, DONE
        end

NbrScale = Trans.NbrScale
npix = Trans.npix
nside = Trans.nside
lmax = Trans.lmax

 
if Trans.MeyerWave EQ 1 then filter=1
if Trans.DifInSH   EQ 1 then filter=1
imag = reform(trans.coef(*,0,*)*0)


if not keyword_set(filter)  then begin 

	for comp=0, 2 do begin 
		Imag(*,comp) = Trans.coef[*,0,comp]
		for j=1, NbrScale-1 do Imag(*,comp) = Imag(*,comp) + Trans.coef[*,j,comp]
	endfor

	mrsp_teb2tqu, Imag, ImagTQU
	Imag = ImagTQU

end else begin

	for comp=0, 2 do begin 
   
		mrs_almtrans, Trans.coef[*,NbrScale-1,comp], ALMTrans, lmax=lmax 
		; help, lmax
		LScale = ALMTrans.ALM
		ech = 2^(NbrScale-2)

		for j=(NbrScale-2), 0, -1 do begin
			if Trans.MeyerWave EQ 0 then begin
				compute_gtilde, lmax, ech, gtilde
				compute_htilde, lmax, ech, htilde
			end else begin
				hgmey, lmax, j, htilde, gtilde, dif=dif
			end

			alm_product2, LScale, htilde, ConvH
			
			mrs_almtrans, Trans.coef[*,j,comp], ALMTrans, lmax=lmax  
      
			WScale = ALMTrans.ALM

			alm_product2, WScale, gtilde, ConvG

			LScale = ConvH + ConvG
			print, "Scale ", j+1, " Re = ", sigma(LScale[*,0]), " IM = ", sigma(LScale[*,1])
			ech = ech/2
		endfor

		ALMTrans.ALM = LScale
		; hs,  ALMTrans
  
		vs = size(ALMTrans.alm)
		n_alm = vs[1]

		if comp EQ 0 then begin
			coef_alm = fltarr( n_alm, 2, 3 )
			coef_alm[*,*,0] = ALMTrans.alm
			ALM = {PixelType: ALMTrans.PixelType, tab: ALMTrans.tab, complex_alm: ALMTrans.complex_alm, nside : ALMTrans.nside, npix: ALMTrans.npix, ALM : coef_alm, norm: ALMTrans.norm, NormVal: ALMTrans.NormVal,lmin: ALMTrans.lmin,lmax: ALMTrans.lmax, TabNbrM: ALMTrans.TabNbrM, index: ALMTrans.index }
		end else ALM.alm[*,*,comp] = ALMTrans.alm
  
		mrsp_almrec, ALM, Imag
   
	endfor ; fin de boucle sur les 3 comp
end

DONE:

END

