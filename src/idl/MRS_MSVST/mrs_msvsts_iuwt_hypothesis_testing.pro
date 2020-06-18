;+
; NAME:
;        mrs_msvsts_iuwt_hypothesis_testing
;
; PURPOSE:
;	Computes the MS-VSTS + Isotropic Undecimate Wavelet Transform of a Poisson Image,
;   perform hypothesis testing on coefficients, returns the multi-resolution support
;   and the denoised image using direct reconstruction
;
;
; CALLING:
;
;   mrs_msvsts_iuwt_hypothesis_testing,image,image_vst,support,image_rec,NbrScale=NbrScale
;   coef_seuil=coef_seuil,First_Scale=First_Scale,background=background
;       
;
; INPUTS:
;     image -- IDL array of healpix map: Input image to be denoised
;    
; OUTPUTS:
;     image_vst -- MS-VST transform of the image
;     support -- multi-resolution support
;     image_rec -- directly reconstructed denoised image
;
; KEYWORDS:
;      NbrScale  : Number of scales (default is 4)
;      coef_seuil  : determines the threshold for the detection of significant coefficients.
;     For each scale i, the threshold is set to coef_seuil*sigma[i] (default is 5)
;      First_Scale  : if > 2, finer wavelet scales are set to 0. (default is 1)
;      background  : if set, substracts a background to the image.
;      
;
; EXTERNAL CALLS:
;       mrs_wttrans
;       mrs_wtrec
;       mrs_msvsts_iuwt_param_computing
;       mrs_msvsts_iuwt_transform
;
;
;
; EXAMPLE:
;
;       Compute the denoising of an image I with default options
;        The result is stored in Output
;               mrs_msvsts_iuwt_hypothesis_testing, Imag, image_vst,support
;         
; HISTORY:
;	Written: Jérémy Schmitt , 2010
;	May, 2010 File creation
;--------------------------------------------------------------------------------------------------------


pro mrs_msvsts_iuwt_hypothesis_testing,image,image_vst,support,image_rec,NbrScale=NbrScale,coef_seuil=coef_seuil,First_Scale=First_Scale,background=background, write=write

if not keyword_set(coef_seuil) then coef_seuil=5
if not keyword_set(coef_pos) then coef_pos=0
if not keyword_set(First_Scale) then First_Scale=1

SizeImage=size(image)
SizeImage=SizeImage[1]

;;;;;;;;;;;;;Calcul des coefficients de la VST;;;;;;;;;;
mrs_msvsts_iuwt_param_computing,NbrScale,c,b,h,tau1,tau2,tau3,sigma

;;;;;;;;;;;;;Transformation en ondelettes et application de la VST sur les coefficients;;;;;;;
mrs_msvsts_iuwt_transform,image,image_vst,NbrScale=NbrScale

if keyword_set(background) then mrs_msvsts_iuwt_transform, background,background_vst,NbrScale=NbrScale

;;;;;;;;;;;;;Tests d'hypothèses et obtention du support de multi-résolution
Support=fltarr(SizeImage,NbrScale-1)
Support[*,*]=0

if coef_seuil LT 5 then begin
 if keyword_set(background) then begin
  ind=where(abs(image_vst.coef[*,0]) - abs(background_vst.coef[*,0]) GT 5*sigma[0],taillesup)
  if taillesup GT 0 then Support[ind,0]=1
 endif else begin ind=where(abs(image_vst.coef[*,0]) GT 5*sigma[0],taillesup)
  if taillesup GT 0 then Support[ind,0]=1
 endelse
endif else begin if keyword_set(background) then begin
  ind=where(abs(image_vst.coef[*,0]) - abs(background_vst.coef[*,0]) GT coef_seuil*sigma[0],taillesup)
  if taillesup GT 0 then Support[ind,0]=1
 endif else begin ind=where(abs(image_vst.coef[*,0]) GT coef_seuil*sigma[0],taillesup)
  if taillesup GT 0 then Support[ind,0]=1
 endelse
endelse

for scale=1,NbrScale-2 do begin
 if keyword_set(background) then begin
  ind=where(abs(image_vst.coef[*,scale]) - abs(background_vst.coef[*,scale]) GT coef_seuil*sigma[scale],taillesup)
  if taillesup GT 0 then Support[ind,scale]=1
 endif else begin ind=where(abs(image_vst.coef[*,scale]) GT coef_seuil*sigma[scale],taillesup)
  if taillesup GT 0 then Support[ind,scale]=1
 endelse
endfor

if (First_Scale GE 2) then begin
 for sc=1,First_Scale-1 do Support[*,sc-1]=0
endif

;reconstruction directe
coef_signif=image_vst.coef
coef_signif[*,0:NbrScale-2]=coef_signif[*,0:NbrScale-2]*Support[*,0:NbrScale-2]

TJaJ=b[NbrScale-1]*sgn(coef_signif[*,NbrScale-1]+c[NbrScale-1])*sqrt(abs(coef_signif[*,NbrScale-1]+c[NbrScale-1]))

tot=total(coef_signif[*,0:NbrScale-2],2)+Tjaj

image_rec=tot^2/(b[0]^2)-c[0]^2

if keyword_set(write) then save,filename='resultat_test_hypothese.xdr',image,image_vst,support,tjaj,tot,image_rec

end