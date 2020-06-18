;+
; NAME:
;        mrs_msvsts_curv_denoising
;
; PURPOSE:
;	Compute Poisson denoising on spherical HEALPix data with MS-VSTS + Curvelets method.
;
;
; CALLING:
;
;    mrs_msvsts_curv_denoising,image,image_reconstruite,NbrScale=NbrScale,niter=niter,HSD=HSD,coef_seuil=coef_seuil,
;    coef_pos=coef_pos,First_Scale=First_Scale,mask=mask,filter=filter,pyr=pyr,background=background,expo=expo
;       
;
; INPUTS:
;     image -- IDL array of healpix map: Input image to be denoised
;     optional : background : if set, substracts a background to the image
;     optional : supentree : if set, use a given multi-resolution support instead of computing it with the procedure ms_vst_tests_hypotheses
;    
; OUTPUTS:
;     image_reconstruite -- IDL array of heapix denoised image
;     optional : support -- multi-resolution support of the image
;     optional : back_reconstruit : if set, returns the reconstructed background (need the keyword separation)
;     
;
; KEYWORDS:
;      NbrScale  : Number of scales (default is 4)
;      niter  : Number of iterations
;      HSD  : if set, the denoised image will be recontructed using the Hybrid Steepest Descent Method (soft thresholding at each 
;     iteration of the reconstruction)
;      coef_seuil  : determines the threshold for the detection of significant coefficients.
;     For each scale i, the threshold is set to coef_seuil*sigma[i] (default is 5)
;      coef_pos  : if set, negative wavelets coefficients are set to 0.
;      First_Scale  : if > 2, finer wavelet scales are set to 0. (default is 1)
;      mask  : if set, enables impainting with the given mask
;      filter  : if set, the inverse wavelet transform will be computed using filters. Else, it will be obtained by a simple addition ;     of all wavelet scales.
;      pyr  : if set, use pyramidal wavelet transform for the soft thresholding
;      expo  : if set, decreases the thresold exponentially at each step of the HSD. Else, decreases the threshold linearly
;      alm : if set, thresholding is made on alm coefficients instead of wavelet coefficients
;      curvelets : if set, thresholding is made on curvelets coefficients instead of wavelet coefficients
;      separation : if set, compute separately the sources and the background
;      update_support : if set, update the multi-resoluation support at each iteration
;      split_support : if set, splits the multi-resolution support
;      
;
; SUBROUTINES
;
;      partie_positive : projection on the non-negative orthant
;      mrs_msvsts_IUWT_reconstruction : iterative reconstruction
;      ms_vst_hypothesis_testing : ms_vsts transform + hypothesis testing
;
;
; EXTERNAL CALLS:
;       ms_vst_test_hypotheses 
;       mrs_wttrans
;       mrs_wtrec
;       mrs_pwttrans
;       mrs_pwtrec
;
; EXAMPLE:
;
;       Compute the denoising of an image I with default options
;        The result is stored in Output
;               ms_vst_recontruction, Imag, Output, NbrScale=6, /hsd
;
;       Compute the denoising of an image I with background extraction and update of multi-resolution_support
;        The result is stored in Output
;               ms_vst_recontruction, Imag, Output, NbrScale=6, /hsd, background=background, /update_support
;
;       Compute the denoising + inpainting of an image I with missing data
;        The result is stored in Output
;               ms_vst_recontruction, Imag, Output, NbrScale=6, /hsd, niter=50, mask=mask
;         
; HISTORY:
;	Written: Jérémy Schmitt & Jean-Luc Starck, 2009
;	February, 2009 File creation
;--------------------------------------------------------------------------------------------------------


pro mrs_msvsts_curv_denoising,image,support,image_rec,scalerec,nbrscale=nbrscale,coef_seuil=coef_seuil

if not keyword_set(coef_seuil) then coef_seuil=5

SizeImage=size(image)
SizeImage=SizeImage[1]

calcul_vst,nbrscale,c,b,h,tau1,tau2,tau3,sigma

mrs_msvsts_curv_transform,image,curvtrans,nbrscale=nbrscale

tabnorm=curvtrans.tabnorm

support=curvtrans

scalerec=fltarr(SizeImage,nbrscale)
scalerec[*,nbrscale-1]=curvtrans.lastscale

for j=1,nbrscale-2 do begin
 nbrscalerid=support.tabnbrscalerid[j]
  for k=0,nbrscalerid-1 do begin
   
   ;my_command = 'Scale'+strcompress(string(j), /remove_all)+'_'+strcompress(string(k), /remove_all)
   ;my_command = my_command + '=mrs_curget(curvtrans,'+strcompress(string(j), /remove_all)+','+strcompress(string(k), /remove_all)+')'
    ;print, 'cmd = ',  my_command
   ;ACK = EXECUTE( my_command) 
   Scalej_k=mrs_curget(curvtrans,j,k)
   
   ;sizescale=size(Scalej_k)
   ;size1=sizescale[1]
   ;size2=sizescale[2]
   ;size3=sizescale[3]
   ;size4=sizescale[4] 
   
   Scalej_k_norma=float(Scalej_k)/tabnorm[j,k]
   
   ind=where(abs(Scalej_k_norma) LE coef_seuil,taillesup)
   if taillesup GT 0 then Scalej_k[ind]=0
   
   mrs_curput,support,Scalej_k,j,k
   
  end

end

nbrscalerid=support.tabnbrscalerid[0]
for k=0,nbrscalerid-1 do begin
   Scale0_k=mrs_curget(curvtrans,0,k)
   
   Scale0_k_norma=float(Scale0_k)/tabnorm[0,k]
   
   ind=where(abs(Scale0_k_norma) LE 10*coef_seuil,taillesup)
   if taillesup GT 0 then Scale0_k[ind]=0
   
   mrs_curput,support,Scale0_k,0,k
end
  


for j=0,nbrscale-2 do begin

 my_command = 'mrs_ridrec,support.ridscale'+strcompress(string(j+1), /remove_all)+',scale'
 ;print,my_command
 ACK = EXECUTE (my_command)
 scalerec[*,j]=scale
 ;scalerec[*,0]=0
 
end

TJaJ=b[NbrScale-1]*sgn(scalerec[*,NbrScale-1]+c[NbrScale-1])*sqrt(abs(scalerec[*,NbrScale-1]+c[NbrScale-1]))

tot=total(scalerec[*,0:NbrScale-2],2)+Tjaj

image_rec=tot^2/(b[0]^2)-c[0]^2

end