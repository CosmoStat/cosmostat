;+
; NAME:
;        mrs_msvsts_curv_denoising
;
; PURPOSE:
;	Compute Poisson denoising on spherical HEALPix data with MS-VSTS + Curvelet Transform method.
;
;
; CALLING:
;
;    mrs_msvsts_curv_denoising,image,image_reconstruite,support,nbrscale=nbrscale,
;    coef_seuil=coef_seuil,suppr_scale1=suppr_scale1,hsd=hsd,niter=niter
;       
;
; INPUTS:
;     image -- IDL array of healpix map: Input image to be denoised
;  
;    
; OUTPUTS:
;     image_reconstruite -- IDL array of heapix denoised image
;     support -- multi-resolution support of the image
;     
;     
;
; KEYWORDS:
;      NbrScale  : Number of scales (default is 4)

;      HSD  : if set, the denoised estimate will be recontructed using the Hybrid Steepest Descent Method
;      (soft thresholding at each iteration of the reconstruction). If not set, the estimate is direclty reconstructed.
;      niter  : Number of iterations for HSD algorithm 
;      coef_seuil  : determines the threshold for the detection of significant coefficients.
;      For each scale i, the threshold is set to coef_seuil*sigma[i] (default is 5)
;      suppr_scale1 : if set, remove the finest scale from the reconstructed estimate
;      
;      
;
; SUBROUTINES
;
;      partie_positive : projection on the non-negative orthant
;
;
; EXTERNAL CALLS:
;       mrs_msvst_curv_transform
;       calcul_vst 
;       mrs_wttrans
;       mrs_wtrec
;       mrs_pwttrans
;       mrs_pwtrec
;
; EXAMPLE:
;
;       Compute the denoising of an image I with direct reconstruction
;        The result is stored in Image_Rec
;               mrs_msvsts_curv_denoising, Image, Image_Rec, Support, NbrScale=6
;
;       Compute the denoising of an image I with iterative reconstruction based on Hybrid Steepest Descent
;        The result is stored in Image_Rec
;               mrs_msvsts_curv_denoising, Image, Image_Rec, Support, NbrScale=6, /hsd
;  
;         
; HISTORY:
;	Written: Jérémy Schmitt & Jean-Luc Starck, 2009
;	February, 2009 File creation
;--------------------------------------------------------------------------------------------------------

function partie_positive,y
taille=size(y)
taille=taille[1]

pos=fltarr(taille)
for j=ulong(0),taille-1 do begin
 if (y[j] LE 0) then pos[j]=0 else pos[j]=y[j]
endfor
return,pos

end

pro mrs_msvsts_curv_denoising,image,image_reconstruite,support,nbrscale=nbrscale,coef_seuil=coef_seuil,suppr_scale1=suppr_scale1,hsd=hsd,niter=niter

if not keyword_set(coef_seuil) then coef_seuil=5

SizeImage=size(image)
SizeImage=SizeImage[1]

mrs_msvsts_IUWT_param_computing,nbrscale,c,b,h,tau1,tau2,tau3,sigma

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
   
   ;Scalej_k_norma=float(Scalej_k)/tabnorm[j,k]
   Scalej_k_norma=float(Scalej_k)/tabnorm[k,j]
   
   ind=where(abs(Scalej_k_norma) LE coef_seuil,taillesup)
   if taillesup GT 0 then Scalej_k[ind]=0
   
   mrs_curput,support,Scalej_k,j,k
   
  end

end

if keyword_set(suppr_scale1) then begin
nbrscalerid=support.tabnbrscalerid[0]
 for k=0,nbrscalerid-1 do begin   
   mrs_curput,support,0,0,k
end
endif else begin nbrscalerid=support.tabnbrscalerid[0]
for k=0,nbrscalerid-1 do begin
   Scale0_k=mrs_curget(curvtrans,0,k)
   
   Scale0_k_norma=float(Scale0_k)/tabnorm[0,k]
   
   ind=where(abs(Scale0_k_norma) LE 10*coef_seuil,taillesup)
   if taillesup GT 0 then Scale0_k[ind]=0
   
   mrs_curput,support,Scale0_k,0,k
end
endelse
  


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

;;reconstruction iterative

if keyword_set(hsd) then begin


x=fltarr(sizeimage,niter+1)
x[*]=0

for i=0.,niter-1 do begin
 print,'n=',i
 residu=image-x[*,i]
 ;tvs,residu+1,/log
 mrs_curtrans,residu,residu_trans,nbrscale=nbrscale,/undec
  for j=0,nbrscale-2 do begin
   nbrscalerid=support.tabnbrscalerid[j]
   for k=0,nbrscalerid-1 do begin
    scjk=mrs_curget(residu_trans,j,k)
    ;print,'avant seuillage'
    ;print,'min(scjk)',min(scjk)
    ;print,'max(scjk)',max(scjk)
    ;print,'mean(scjk)',mean(scjk)
    ;print,'sigma(scjk)',sigma(scjk)
    ;print,'total(scjk)',total(scjk)
    supjk=mrs_curget(support,j,k)
    ;ind=where(abs(supjk) GT 0,taillesup)
    ;if taillesup GT 0 then supjk[ind]=1.
    ;scjk=scjk*supjk
    ind=where(abs(supjk) EQ 0,taillesup)
    if taillesup GT 0 then scjk[ind]=0.
    ;print,'apres seuillage'
    ;print,'min(scjk)',min(scjk)
    ;print,'max(scjk)',max(scjk)
    ;print,'mean(scjk)',mean(scjk)
    ;print,'sigma(scjk)',sigma(scjk)
    ;print,'total(scjk)',total(scjk)
    mrs_curput,residu_trans,scjk,j,k
   end
  end
 mrs_currec,residu_trans,residu_rec
 ;tvs,residu_rec+1,/log
 x[*,i+1]=x[*,i]+residu_rec
 
 seuil=1.*(niter-(i+1.))/(niter-1.)
 print, 'seuil=',seuil
 mrs_wttrans,x[*,i+1],w,nbrscale=nbrscale
 for l=0,nbrscale-1 do begin
  sig=sigma(w.coef[*,l])
  thresh=sig*seuil
  coeff=w.coef[*,l]
  softthreshold,coeff,thresh
  w.coef[*,l]=coeff
 endfor
 mrs_wtrec,w,wrec
 x[*,i+1]=wrec
 x[*,i+1]=partie_positive(x[*,i+1])
 save,filename='resultat_reconstruction_curv.xdr',image,x,image_vst,support,NbrScale,niter
endfor
image_reconstruite=x[*,niter]

endif else image_reconstruite=image_rec

end