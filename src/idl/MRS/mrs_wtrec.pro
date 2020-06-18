;+
; NAME:
;        mrs_wtrec
;
; PURPOSE:
;	Compute the inverse wavelet transform on the sphere.
;
;
; CALLING:
;
;     mrs_wtrec, Trans, Rec, filter=filter        
;    
; INPUT:
;     Trans -- IDL structures with the following fields:  
;                  NbrScale : int = number of scales 
;                     nside : int = Healpix nside parameter (0 for a Glesp image)
;                      lmax : int = maximum l value in the Spherical Harmonic Space (Healpix)
;                      npix : long = Number of pixels of the input image (12*nside*nside)
;                      Coef : fltarr[npix,NbrScale] = wavelet transform of the data
;                             Coef[*,0] = wavelet coefficients of the finest scale (highest frequencies).
;                             Coef[*,NbrScale-1] = coarsest scale (lowest frequencies). 
;                      lmax : int = lmax parameter at the first scale
;		 Healpix_with_Glesp : int = 1 if the keyword Healpix_with_Glesp used, otherwise 0
;				   UseGLESP : int = 1 if the input image was in Glesp format, otherwise 0
;				  MeyerWave : int = 1 if the keyword MeyerWave used, otherwise 0
;					DifInSH : int = 1 if the keyword DifInSH used, otherwise 0
;						 nx : int = number of rings (Glesp parameter), otherwise, 0
;						 np : int = max number of pixel on a ring (Glesp parameter), otherwise, 0
;					  x_sky : fltarray with COS( THETA ) for each ring (Glesp parameter), otherwise, 0
;					  y_sky : long 1D array number of pixels/ring (Glesp parameter), otherwise, 0
;
; OUTPUT:
;		Imag -- IDL array of healpix map: reconstructed image from the wavelet coefficients or Glesp image IDL structure if UseGlesp=1
;
; KEYWORDS:
;      filter : Use filters for the reconstructions. If this keyword is not set, the reconstructed image
;               is obtained by a simple addition of all wavelet scales. Automaticaly applied if keyword 
;				MeyerWave, NeedletWave or DifInSH were set at the wavelet decomposition.
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
;               mrs_wtrec, Trans, Imag 
;         
; HISTORY:
;	Written: Pierrick Abrial & Jean-Luc Starck, 2004
;	December, 2004 File creation
;------------------------------------------------------------------------------

pro mrs_wtrec, Trans, Imag, filter=filter 

 
if N_PARAMS() LT 2 or N_PARAMS() GE 3 then begin 
        print, 'CALLING SEQUENCE: mrs_wtrec, Trans, Imag, filter=filter'
        goto, DONE
        end

NbrScale = Trans.NbrScale
npix =  Trans.npix
nside = Trans.nside
lmax = Trans.lmax
UseGlesp = Trans.USEGLESP
if UseGlesp EQ  1 then Gmap = {T_SKY: Trans.coef[*,0], x_sky:Trans.x_sky, y_sky:Trans.y_sky, nx:Trans.nx, np:Trans.np}
 
if Trans.MeyerWave EQ 1 then filter=1
if Trans.DifInSH   EQ 1 then filter=1
if Trans.NeedletWave EQ 1 then filter=1

; Simple summation for the reconstruction
if not keyword_set(filter) or Trans.Healpix_with_Glesp EQ 1 then begin 
   ; Imag = Trans.coef[*,0]
   ; for j=1,NbrScale-1 do Imag = Imag + Trans.coef[*,j]
   Imag = total(double(Trans.coef),2)
   if UseGlesp EQ  1 then begin
      Gmap.T_SKY = Imag
      Imag = Gmap
   end
end else begin
; Reconstruction in the spherical harmonic domain
  if UseGlesp EQ  0 then mrs_almtrans, Trans.coef[*,NbrScale-1], ALMTrans, lmax=lmax $
  else begin
     Gmap.T_SKY = Trans.coef[*,NbrScale-1]
     mrs_almtrans, Gmap, ALMTrans, lmax=lmax 
  end
  ; help, lmax
  LScale = ALMTrans.ALM
  ech = 2^(NbrScale-2)
  for j =(NbrScale-2),0,-1 do begin
      if Trans.MeyerWave EQ 1  then begin ; Meyer filters
                hgmey, lmax, j, htilde, gtilde, dif=dif
      end else  if Trans.NeedletWave EQ 1 then  begin ; Needet filters
              gtilde = Trans.tabpsi[*,j]
              htilde = Trans.tabphi[*,NbrScale-2]
              if j NE NbrScale-2 then htilde[*] = 1
     end else begin  ; B3-spline filters
        compute_gtilde, lmax, ech, gtilde
        compute_htilde, lmax, ech, htilde
      end       
        
      alm_product2, LScale, htilde, ConvH
      
      if UseGlesp EQ  0 then mrs_almtrans, Trans.coef[*,j], ALMTrans, lmax=lmax  $
      else begin
         Gmap.T_SKY = Trans.coef[*,j]
         mrs_almtrans, Gmap, ALM, lmax=lmax 
      end
      WScale = ALMTrans.ALM
      alm_product2, WScale, gtilde, ConvG
      LScale = ConvH + ConvG
      ; print, "Scale ", j+1, " Re = ", sigma(LScale[*,0]), " IM = ", sigma(LScale[*,1])
      ech = ech/2
   endfor
   ALMTrans.ALM = LScale
   ; hs,  ALMTrans
   mrs_almrec, ALMTrans, Imag  
end

DONE:

END

