;+
; NAME:
;        mrs_msvsts_multichannel_deconvolution
;
; PURPOSE:
;	Compute multichannel Poisson deconvolution on spherical 2D-1D HEALPix data with MS-VSTS + multichannel Wavelet Transform method.
;
;
; CALLING:
;
;    mrs_msvsts_multichannel_denoising,input,solution,NbrScale1=NbrScale1, NbrScale2=NbrScale2,niter=niter,beam=beam,regularization=regularization
;       
;
; INPUTS:
;     input -- (IDL array) multichannel healpix map: Input data to be denoised
;     beam -- (IDL array) multichannel convolution beam
;       
; OUTPUTS:
;     Solution -- (IDL array) multichannel HEALPix denoised data
;     
;
; KEYWORDS:
;      NbrScale1 -- Number of scales for the two spatial dimensions (default is 6) 
;      NbrScale2 -- Number of scales for the non-spatial dimension (time or energy) (default is 6) 
;      niter -- Number of iterations 
;      regularization -- Use of a regularization parameter. Improves the speed of the algorithm

;      
;
; EXAMPLE:
;
;       Compute the deconvolution of a multichannel data with a given
;       multichannel beam
;        The result is stored in Output
;               mrs_msvsts_multichannel_denoising, Data, Output, NbrScale1=6, NbrScale2=6, niter=100, beam=beam
;
;       Compute the deconvolution of a multichannel data with a given
;       multichannel beam and a regularization parameter
;        The result is stored in Output
;               mrs_msvsts_multichannel_denoising, Data, Output, NbrScale1=6, NbrScale2=6, niter=100, beam=beam, regularization=regularization

;         
; HISTORY:
;	Written: JŽrŽmy Schmitt
;	November, 2011 File creation
;--------------------------------------------------------------------------------------------------------


pro calcul_sigma_multichannel,sizemc,nbrscale1,nbrscale2,sigmadd,sigmadc,sigmacd

test=fltarr(196608,sizemc)
test[*]=100
testpoisson=poisson_image(test)

msvsts_multichannel,testpoisson,solution,trans,support,image_wt_mc,scaling_wt_mc,cdxztab,dcxztab,ddxztab,NbrScale1=NbrScale1,NbrScale2=NbrScale2,niter=1

sigmadd=fltarr(sizemc,nbrscale1,nbrscale2)
sigmadc=fltarr(sizemc,nbrscale1)
sigmacd=fltarr(sizemc,nbrscale2)

for mc=0,sizemc-1 do begin
 for sx=0,nbrscale1-1 do begin
  for sz=0,nbrscale2-1 do begin
   sigmadd[mc,sx,sz]=sigma(ddxztab[*,mc,sx,sz])
   sigmadc[mc,sx]=sigma(dcxztab[*,mc,sx])
   sigmacd[mc,sz]=sigma(cdxztab[*,mc,sz])
  end
 end
end

save,filename='resultat_calcul_sigma_multichannel.xdr',test,solution,sigmadd,sigmadc,sigmacd,testpoisson,cdxztab,dcxztab,ddxztab

end

pro generer_donnees_mc,data,sizeMC
 chargement_donnees,model,sources,modelsources,imagepoisson
 sizeimage=size(imagepoisson)
 sizeimage=sizeimage[1]
 data=fltarr(sizeimage,sizemc)
 for j=0.,sizemc-1 do begin
  modelj=modelsources/(j+1.)
  imagej=poisson_image(modelj)
  data[*,j]=imagej
 end
end

function devlpl,a,n,x
 term = a[n-1]
 for i=0,n-1-1 do begin
  j=n-1-1-i
  term = a[j]+term*x
 end
 return,term
end

function stvaln,p
 xden = [0.993484626060e-1,0.588581570495e0,0.531103462366e0,0.103537752850e0,0.38560700634e-2]
 xnum = [-0.322232431088e0,-1.000000000000e0,-0.342242088547e0,-0.204231210245e-1,-0.453642210148e-4]
 K1=5
 if (p GT 0.5) then begin
  sign=1.
  z=1.-p
 endif else begin
  sign=-1.
  z=p
 endelse
 y=sqrt(-(2.*alog(z)))
 valn=y+devlpl(xnum,K1,y)/devlpl(xden,K1,y)
 valn=sign*valn
 return,valn
end

function criticalThreshGauss,tailProba
 tab1=[1.-tailProba,0.5]
 max1=max(tab1)
 tab2=[1.,max1]
 min1=min(tab2)
 valn=stvaln(min1)
 tab3=[0,valn]
 max2=max(tab3)
 return,max2
end

pro gaussHardThresholdDenoise,waveletData,alpha,sigma,support,j,k
 SizeDataMC=size(waveletData)
 SizeImage=SizeDataMC[1]
 SizeMC=SizeDataMC[2]
 ;cthresh=criticalThreshGauss(alpha/2.)
 ;cthresh=4.
 cthresh=6.
 if (j EQ 0) then cthresh=10.
 thresh=cthresh*sigma
 ;essai
 ;if j EQ 0 then thresh=50.*thresh
 ;ind=where(abs(waveletdata[*,*,j,k]) GT thresh,taillesup)
 
 
 
 ;ind=where(abs(waveletdata) GT thresh,taillesup)
 
 ;;;;;coeffs positifs
 ;ind=where(waveletdata GT thresh,taillesup)
 
 ;3 premieres echelles : coef pos
 ;if (j GE 2) then ind=where(abs(waveletdata) GT thresh,taillesup) else ind=where(waveletdata GT thresh,taillesup)
 ind=where(abs(waveletdata) GT thresh,taillesup)
 
 suppjk=fltarr(sizeimage,sizemc)
 if taillesup GT 0 then begin
  suppjk[ind]=1
  Support[*,*,j,k]=suppjk[*,*]
 end
 
 ;;temporaire
   ;print,'j=',j,'k=',k
   ;print,'info ind'
   ;info,ind
   ;print,'info coeff'
   ;info,abs(waveletdata)
   ;print,'thresh = ',thresh

end

function partie_positive,y
 pos=y
 ind=where(y LT 0,taillesup)
 if (taillesup GT 0) then pos[ind]=0
 return,pos
end


pro calcul_variancez,scalex,scalez,variance
 tau1=[1.00000,1.00003,1.00008,1.00018,1.00024,1.00037,1.00067,1.00063,1.00032,1.00001,1.00304]
 tau2=[1.00000, 0.0225476, 0.00462573, 0.00111173, 0.000278866, 7.16454e-05, 1.90052e-05, 5.34395e-06, 1.67009e-06, 6.20837e-07, 3.17044e-07]
 tau3=[1.00000, 0.000637954, 2.74081e-05, 1.58719e-06, 9.95906e-08, 6.52659e-09, 4.52801e-10, 3.48761e-11, 3.26746e-12, 4.31919e-13, 1.01216e-13]
 t1=[1,1,1,1,1,1,1,1,1,1,1]
 t2=[1, 0.2734375, 0.12347412109375, 0.0603594779968262, 0.0300147198140621, 0.014986945054261, 0.00749092722048772, 0.00374514565088013, 0.00187253308688793, 0.00093626157632392, 0.000468130167278172]
 t3=[1, 0.08447265625, 0.017338752746582, 0.00414866860955954, 0.00102616434651281, 0.00025586256661736, 6.39233749606246e-5, 1.59782042631771e-5, 3.99438613268932e-6, 9.9858622538761e-7, 2.49645912118706e-7]
 sumhh=[0, 0.375, 0.15771484375, 0.0760841369628906, 0.0377178490161896, 0.0188190417829901, 0.00940455525596917, 0.00470165753587537, 0.00235075127553241, 0.00117536595181247, 0.000587681765180673, 0.000293840731250224]
 ;tt21=tau2[scalex]^2*t2[scalez-1]
 ;tt22=tau2[scalex]^2*t2[scalez]
 ;shh=tau2[scalex]^2*sumhh[scalez]
 tt21=tau2[scalex]*t2[scalez-1]
 tt22=tau2[scalex]*t2[scalez]
 shh=tau2[scalex]*sumhh[scalez]
 variance=tt21/4.+tt22/4.-shh/2.
end

pro calcul_variancex,scalex,scalez,variance
 tau1=[1.00000,1.00003,1.00008,1.00018,1.00024,1.00037,1.00067,1.00063,1.00032,1.00001,1.00304]
 tau2=[1.00000, 0.0225476, 0.00462573, 0.00111173, 0.000278866, 7.16454e-05, 1.90052e-05, 5.34395e-06, 1.67009e-06, 6.20837e-07, 3.17044e-07]
 tau3=[1.00000, 0.000637954, 2.74081e-05, 1.58719e-06, 9.95906e-08, 6.52659e-09, 4.52801e-10, 3.48761e-11, 3.26746e-12, 4.31919e-13, 1.01216e-13]
 t1=[1,1,1,1,1,1,1,1,1,1,1]
 t2=[1, 0.2734375, 0.12347412109375, 0.0603594779968262, 0.0300147198140621, 0.014986945054261, 0.00749092722048772, 0.00374514565088013, 0.00187253308688793, 0.00093626157632392, 0.000468130167278172]
 t3=[1, 0.08447265625, 0.017338752746582, 0.00414866860955954, 0.00102616434651281, 0.00025586256661736, 6.39233749606246e-5, 1.59782042631771e-5, 3.99438613268932e-6, 9.9858622538761e-7, 2.49645912118706e-7]
 sumhh=[0, 0.375, 0.15771484375, 0.0760841369628906, 0.0377178490161896, 0.0188190417829901, 0.00940455525596917, 0.00470165753587537, 0.00235075127553241, 0.00117536595181247, 0.000587681765180673, 0.000293840731250224]
 ;tt21=tau2[scalex-1]^2*t2[scalez]
 ;tt22=tau2[scalex]^2*t2[scalez]
 tt21=tau2[scalex-1]*t2[scalez]
 tt22=tau2[scalex]*t2[scalez]
 ;shh=t2[scalex]^2*t2[scalez]
 ;shh=sumhh[scalex]^2*t2[scalez]
 shh=sumhh[scalex]*t2[scalez]
 variance=tt21/4.+tt22/4.-shh/2.
end



pro convolz,input,output,step

B3=[1./16., 1./4., 3./8., 1./4., 1./16.]
h=B3

filterLen=5
SizeDataMC=size(input)
SizeImage=SizeDataMC[1]
SizeMC=SizeDataMC[2]
output=fltarr(SizeImage,SizeMC)

;for x=ulong(0),sizeImage-1 do begin
 for z=0,sizeMC-1 do begin
  ;output[x,z]=0
  output[*,z]=0
  for p=0,filterLen-1 do begin
   c=(filterLen/2-p)*step
   index=z-c
   ;inp=input
   
   ;;get_index_mirror
   ;if (index LT 0) then begin
   ; index=-index
   ; if (index GE SizeMC) then index=SizeMC-1
   ;endif else begin
   ; if (index GE SizeMC) then index=2*(SizeMC-1)-index
   ; if (index LT 0) then index=0
   ;endelse
   
   ;get_index_period
   ;if (index LT 0) then begin
   ; while (index LT 0) do index = index+SizeMC
   ;endif else begin
   ; if (index GE SizeMC) then begin 
   ;  while (index GE SizeMC) do index=index-SizeMC
   ; endif
   ;endelse
   
   ;get_index_cont
   ;if (index LT 0) then index=0
   ;if (index GE SizeMC) then index=sizemc-1
   
   ;get_index_zero
   if (index GE 0) AND (index LT SizeMC) then begin
    output[*,z]=output[*,z]+input[*,index]*h[filterLen-p-1.]
   end
   ;output[x,z]=output[x,z]+input[x,index]*h[filterLen-p-1.]
   ;output[*,z]=output[*,z]+input[*,index]*h[filterLen-p-1.]
   
  end
 end
;end

end

pro MSVST,data,vstdata,scalex,scalez

SizeDataMC=size(data)
SizeImage=SizeDataMC[1]
SizeMC=SizeDataMC[2]
vstdata=fltarr(SizeImage,SizeMC)

tau1=[[1.00001, 1.00004, 1.00012, 1.00002, 1.00001, 1.00000, 0.999994, 1.00019],$
      [1.62503, 1.62507, 1.62520, 1.62502, 1.62494, 1.62501, 1.62499, 1.62506],$      
      [0.687508, 0.687524, 0.687579, 0.687512, 0.687506, 0.687502, 0.687496, 0.687561],$      
      [1.06251, 1.06254, 1.06262, 1.06252, 1.06251, 1.06250, 1.06249, 1.06256],$     
      [1.06252, 1.06254, 1.06261, 1.06252, 1.06250, 1.06250, 1.06249, 1.06254],$     
      [1.68753, 1.68758, 1.68771, 1.68752, 1.68745, 1.68751, 1.68749, 1.68762],$     
      [2.93756, 2.93764, 2.93786, 2.93754, 2.93738, 2.93751, 2.93748, 2.93761],$      
      [3.43757, 3.43766, 3.43792, 3.43754, 3.43736, 3.43751, 3.43747, 3.43761]]
      
tau2=[[0.0227758,   0.00472934,   0.00116309,  0.000305187,  8.55772e-05,  2.67357e-05,  9.93561e-06,  5.08527e-06],$      
      [0.0201068,   0.00417513,   0.00102679,  0.000269423,  7.55486e-05,  2.35995e-05,  8.77128e-06,  4.48985e-06],$     
      [0.00471531,  0.000979122,  0.000240797,  6.31833e-05,  1.77172e-05,  5.53520e-06,  2.05698e-06,  1.05295e-06],$     
      [0.00685055,   0.00142250,  0.000349837,  9.17946e-05,  2.57401e-05,  8.04173e-06,  2.98844e-06,  1.52983e-06],$      
      [0.00524912,   0.00108997,  0.000268057,  7.03362e-05,  1.97229e-05,  6.16191e-06,  2.28984e-06,  1.17229e-06],$     
      [0.00951959,   0.00197672,  0.000486137,  0.000127559,  3.57686e-05,  1.11733e-05,  4.15277e-06,  2.12571e-06],$      
      [0.0201958,   0.00419360,   0.00103134,  0.000270615,  7.58829e-05,  2.37040e-05,  8.81009e-06,  4.50972e-06],$      
      [0.0244662,   0.00508035,   0.00124942,  0.000327838,  9.19286e-05,  2.87162e-05,  1.06730e-05,  5.46332e-06]] 
      
tau3=[[0.000649356,  2.84927e-05,  1.71803e-06,  1.16725e-07,  8.94354e-09,  8.36870e-10,  1.10589e-10,  2.58751e-11],$ 
  [0.000337044,  1.47890e-05,  8.91731e-07,  6.05851e-08,  4.64209e-09,  4.34371e-10,  5.74004e-11,  1.34307e-11],$ 
  [4.45481e-05,  1.95470e-06,  1.17863e-07,  8.00771e-09,  6.13557e-10,  5.74122e-11,  7.58678e-12,  1.77486e-12],$ 
  [6.92794e-05,  3.03988e-06,  1.83295e-07,  1.24533e-08,  9.54179e-10,  8.92852e-11,  1.17986e-11,  2.76007e-12],$ 
  [5.59625e-05,  2.45555e-06,  1.48062e-07,  1.00595e-08,  7.70770e-10,  7.21226e-11,  9.53072e-12,  2.22987e-12],$ 
  [8.03768e-05,  3.52681e-06,  2.12656e-07,  1.44481e-08,  1.10703e-09,  1.03587e-10,  1.36886e-11,  3.20290e-12],$ 
  [0.000178668,  7.83968e-06,  4.72710e-07,  3.21164e-08,  2.46079e-09,  2.30262e-10,  3.04282e-11,  7.11966e-12],$ 
  [0.000217984,  9.56483e-06,  5.76731e-07,  3.91837e-08,  3.00229e-09,  2.80932e-10,  3.71240e-11,  8.68637e-12]]      


;tau1=[1.00000,1.00003,1.00008,1.00018,1.00024,1.00037,1.00067,1.00063,1.00032,1.00001,1.00304]
;tau2=[1.00000, 0.0225476, 0.00462573, 0.00111173, 0.000278866, 7.16454e-05, 1.90052e-05, 5.34395e-06, 1.67009e-06, 6.20837e-07, 3.17044e-07]
;tau3=[1.00000, 0.000637954, 2.74081e-05, 1.58719e-06, 9.95906e-08, 6.52659e-09, 4.52801e-10, 3.48761e-11, 3.26746e-12, 4.31919e-13, 1.01216e-13]


;t1=[1,1,1,1,1,1,1,1,1,1,1]
;t2=[1, 0.2734375, 0.12347412109375, 0.0603594779968262, 0.0300147198140621, 0.014986945054261, 0.00749092722048772, 0.00374514565088013, 0.00187253308688793, 0.00093626157632392, 0.000468130167278172]
;t3=[1, 0.08447265625, 0.017338752746582, 0.00414866860955954, 0.00102616434651281, 0.00025586256661736, 6.39233749606246e-5, 1.59782042631771e-5, 3.99438613268932e-6, 9.9858622538761e-7, 2.49645912118706e-7]
;tt1=tau1[scalex]^2*t1[scalez]
;tt2=tau2[scalex]^2*t2[scalez]
;tt3=tau3[scalex]^2*t3[scalez]
;tt1=tau1[scalex]*t1[scalez]
;tt2=tau2[scalex]*t2[scalez]
;tt3=tau3[scalex]*t3[scalez]
tt1=tau1[scalex,scalez]
tt2=tau2[scalex,scalez]
tt3=tau3[scalex,scalez]

c=7.*tt2/(8.*tt1)-tt3/(2.*tt2)
b=2.*sqrt(tt1/tt2)


;for j=ulong(0),SizeImage*SizeMC-1 do vstdata[j]=b*sqrt(data[j]+c)
for j=ulong(0),SizeImage-1 do begin
 ;for k=0,SizeMC-1 do vstdata[j,k]=b*sqrt(data[j,k]+c)
 for k=0,SizeMC-1 do vstdata[j,k]=sqrt(data[j,k]+c)
end

end

pro transformz,data,appr,detail,scalex,scalez
 step=(scalez-1.)^2
 ;step=(scalez)^2
 ;step=1
 convolz,data,appr,step
 ;d1=data
 ;for st=0,step do begin
 ; convolz,d1,d2,st
 ; d1=d2
 ;end
 ;appr=d2
 
 detail=data-appr
 ;convolz,appr,detail1,step
 ;d1=appr
 ;for st=0,step do begin
 ; convolz,d1,d2,st
 ; d1=d2
 ;end
 ;detail1=d2
 
 ;detail=data-detail1 
end

pro reconsz,appr,detail,data,scalex,scaley
 ;step=(scalez-1)^2
 ;convolz,appr,data,step
 ;data=data+detail
 data=appr+detail
end

pro transformz_d,approxpxpz,approxpxcz,approxcxpz,approxcxcz,da,dd,scalex,scalez
 step=(scalez-1.)^2
 ;step=(scalez)^2
 ;step=1
 da=approxpxcz-approxcxcz
 dd=approxpxpz-approxcxpz
 convolz,dd,temp,step
; d1=dd
; for st=0,step do begin
;  convolz,d1,d2,st
;  d1=d2
; end
; temp=d2
 dd=dd-temp
end

pro transformz_approx,data,appr,vstdetail,scalex,scalez

SizeDataMC=size(data)
SizeImage=SizeDataMC[1]
SizeMC=SizeDataMC[2]
step=(scalez-1.)^2
;step=(scalez)^2

convolz,data,appr,step
;d1=data
;for st=0,step do begin
; convolz,d1,d2,st
; d1=d2
;end
;appr=d2
MSVST,data,vstdetail,scalex,scalez-1
;MSVST,data,vstdetail,scalex,scalez
MSVST,appr,temp,scalex,scalez
;for j=ulong(0),SizeImage*SizeMC-1 do vstdetail[j]=vstdetail[j]-temp[j]
vstdetail=vstdetail-temp

end

pro transformz_detail,approxpxpz,approxpxcz,approxcxpz,approxcxcz,vstda,vstdd,scalex,scalez

step=(scalez-1.)^2
;step=(scalez)^2
;step=1
MSVST,approxpxcz,vst1,scalex-1,scalez
MSVST,approxcxcz,vst2,scalex,scalez
;MSVST,approxpxcz,vst1,scalex,scalez
;MSVST,approxcxcz,vst2,scalex,scalez
vstda=vst1-vst2

MSVST,approxpxpz,vst1,scalex-1,scalez-1
MSVST,approxcxpz,vst2,scalex,scalez-1
;MSVST,approxpxpz,vst1,scalex,scalez
;MSVST,approxcxpz,vst2,scalex,scalez
vstdd=vst1-vst2
convolz,vstdd,vst1,step
;d1=vstdd
;for st=0,step do begin
; convolz,d1,d2,st
; d1=d2
;end
;vst1=d2
vstdd=vstdd-vst1

end

pro waveletz_approx,data,appr,vstdetail,scalex,scalez

SizeDataMC=size(data)
SizeImage=SizeDataMC[1]
SizeMC=SizeDataMC[2]
step=(scalez-1.)^2
;step=(scalez)^2

convolz,data,appr,step
;d1=data
;for st=0,step do begin
; convolz,d1,d2,st
; d1=d2
;end
;appr=d2
;MSVST,data,vstdetail,scalex,scalez-1
;MSVST,data,vstdetail,scalex,scalez
;MSVST,appr,temp,scalex,scalez
;for j=ulong(0),SizeImage*SizeMC-1 do vstdetail[j]=vstdetail[j]-temp[j]
vstdetail=data-appr

end


pro waveletz_detail,approxpxpz,approxpxcz,approxcxpz,approxcxcz,da,dd,scalex,scalez
 step=(scalez-1.)^2
 ;step=(scalez)^2
 ;step=1
 da=approxpxcz-approxcxcz
 dd=approxpxpz-approxcxpz
 convolz,dd,temp,step
; d1=dd
; for st=0,step do begin
;  convolz,d1,d2,st
;  d1=d2
; end
; temp=d2
 dd=dd-temp
end

pro estimation_support,data,support,nbrscale1=nbrscale1,nbrscale2=nbrscale2

 if not keyword_set(nbrscale1) then nbrscale1=6
 if not keyword_set(nbrscale2) then nbrscale2=6

 ;;B3-Spline
 B3=[1./16., 1./4., 3./8., 1./4., 1./16.]
 h=B3
 
 ;;Default cut p-value at 3.5*sigma
 alpha=0.000465
 
 ;;values of sigma (DD, DC and CD coefficients)
 sigmadd=[[ 1.08061e-09,      0.00000,      0.00000,      0.00000,      0.00000,      0.00000],$
     [0.349935,    0.0398326,    0.0171225,   0.00826525,   0.00410096,   0.00196823],$
     [0.178645,    0.0204763,   0.00878452,   0.00422461,   0.00210050,   0.00110088],$
    [0.0801816,   0.00916309,   0.00393344,   0.00191854,  0.000959409,  0.000484630],$
    [0.0534391,   0.00609373,   0.00257950,   0.00127993,  0.000651180,  0.000358986],$
    [0.0369379,   0.00421229,   0.00177550,  0.000887714,  0.000452417,  0.000258144]]
sigmadc=[0.0361890,   0.00412712,   0.00173962,  0.000869778,  0.000443276,  0.000252928]
sigmacd=[2.61089e-07,   0.00166342,  0.000794248,  0.000437837,  0.000349688,  0.000313635]

 SizeDataMC=size(data)
 SizeImage=SizeDataMC[1]
 SizeMC=SizeDataMC[2]
 image_wt_mc=fltarr(sizeimage,nbrscale1+1,sizemc)
 scaling_wt_mc=fltarr(sizeimage,nbrscale1+1,sizemc)
 vecMC=fltarr(SizeMC)
 trans=fltarr(sizeimage,nbrscale1,sizemc)

 ;coeff_ApAp=fltarr(SizeImage,SizeMC,nbrscale1,nbrscale2)
 ;coeff_AA=fltarr(SizeImage,SizeMC,nbrscale1,nbrscale2)
 ;coeff_AD=fltarr(SizeImage,SizeMC,nbrscale1,nbrscale2)
 ;coeff_DA=fltarr(SizeImage,SizeMC,nbrscale1,nbrscale2)
 ;coeff_DD=fltarr(SizeImage,SizeMC,nbrscale1,nbrscale2)

 support=fltarr(SizeImage,SizeMC,nbrscale1+1,nbrscale2+1)
 support[*]=0

 
 
 cx=fltarr(SizeImage,SizeMC)
 dx=fltarr(SizeImage,SizeMC)
 ccxz=fltarr(SizeImage,SizeMC,nbrscale2+1)
 ccpxz=fltarr(SizeImage,SizeMC,nbrscale2+1)
 cdxz=fltarr(SizeImage,SizeMC)
 dcxz=fltarr(SizeImage,SizeMC)
 ddxz=fltarr(SizeImage,SizeMC)
 temp=fltarr(SizeImage,SizeMC)
 
 cdxztab=fltarr(SizeImage,SizeMC,nbrscale2)
 dcxztab=fltarr(SizeImage,SizeMC,nbrscale1)
 ddxztab=fltarr(SizeImage,SizeMC,nbrscale1,nbrscale2)
 cdxztabp=fltarr(SizeImage,SizeMC,nbrscale2)
 dcxztabp=fltarr(SizeImage,SizeMC,nbrscale1)
 ddxztabp=fltarr(SizeImage,SizeMC,nbrscale1,nbrscale2)

 for j=0,SizeMC-1 do begin
  image=data[*,j]
  mrs_wttrans,image,image_wt,NbrScale=NbrScale1+1
  image_wt_mc[*,*,j]=image_wt.coef
 end
 
 scaling_wt_mc=image_wt_mc
 
 ;for sx=0,nbrscale1-2 do begin
 for sx=0,nbrscale1-1 do begin
  ;;;calcul des coefficients d'échelle de la 1e transformée
  scaling_wt_mc[*,sx,*]=total(image_wt_mc[*,sx:nbrscale1,*],2)
 end

 temp=data
 ccpxz[*,*,0]=temp
 for sz=1,nbrscale2 do begin 
  transformz_approx,temp,cp,cdxyz,0,sz 
  ccpxz[*,*,sz]=cp
  temp=cp
 end
 
 for j=1,nbrscale1 do begin
  print,'scalex =',j
  cx[*,*]=scaling_wt_mc[*,j,*]
  temp=cx
  ccxz[*,*,0]=temp
  

  for k=1,nbrscale2 do begin
   print,'scalez =',k
   ;;calcul des coefficients approx_approx et approx_detail
   transformz_approx,temp,ccx,cdxz,j,k
   ccxz[*,*,k]=ccx
   temp=ccxz[*,*,k]
   ;;Denoising sur les coeff_AD pour j=nbrscale1
   if (j EQ nbrscale1) then begin
    ;calcul_variancez,j,k,variance
    ;calcul_variancez,j,1,variance
    ;sigma=sqrt(variance)
    sigma=sigmacd[k-1]
    gaussHardThresholdDenoise,cdxz,alpha,sigma,support,j,k-1
    cdxztab[*,*,k-1]=cdxz
   end
  end
  
  for k=1,nbrscale2 do begin
   ;;calcul des coefficients detail_approx et detail_detail
   ccpxzprev=fltarr(SizeImage,SizeMC)
   ccpxzact=fltarr(SizeImage,SizeMC)
   ccxzprev=fltarr(SizeImage,SizeMC)
   ccxzact=fltarr(SizeImage,SizeMC)
   ccpxzprev[*,*]=ccpxz[*,*,k-1]
   ccpxzact[*,*]=ccpxz[*,*,k]
   ccxzprev[*,*]=ccxz[*,*,k-1]
   ccxzact[*,*]=ccxz[*,*,k]
   transformz_detail,ccpxzprev,ccpxzact,ccxzprev,ccxzact,dcxz,ddxz,j,k
   ;;Denoising sur les coeff_DA pour k=nbrscale2-1
   if (k EQ nbrscale2) then begin
    ;calcul_variancex,j,k,variance
    ;calcul_variancex,j,0,variance
    ;sigma=sqrt(variance)
    sigma=sigmadc[j-1]
    gaussHardThresholdDenoise,dcxz,alpha,sigma,support,j-1,k
    dcxztab[*,*,j-1]=dcxz
    end
   ;;Denoising sur les coeff_DD
   sigma=sigmadd[j-1,k-1]
   gaussHardThresholdDenoise,ddxz,alpha,sigma,support,j-1,k-1
   ddxztab[*,*,j-1,k-1]=ddxz
  end
  
 for k=0,nbrscale2 do ccpxz[*,*,k]=ccxz[*,*,k]
 
 end


end

pro mrs_msvsts_multichannel_deconvolution,donnees,solution,trans,support,image_wt_mc,scaling_wt_mc,cdxztab,dcxztab,ddxztab,NbrScale1=NbrScale1,NbrScale2=NbrScale2,niter=niter,beam=beam,regularization=regularization
  
 ;keyword beam for multi-scale deconvolution

 if not keyword_set(nbrscale1) then nbrscale1=6
 if not keyword_set(nbrscale2) then nbrscale2=6

 ;;B3-Spline
 B3=[1./16., 1./4., 3./8., 1./4., 1./16.]
 h=B3
 
 ;;Default cut p-value at 3.5*sigma
 alpha=0.000465
 
 ;;values of sigma (DD, DC and CD coefficients)
sigmadd=[[ 1.08061e-09,      0.00000,      0.00000,      0.00000,      0.00000,      0.00000],$
     [0.349935,    0.0398326,    0.0171225,   0.00826525,   0.00410096,   0.00196823],$
     [0.178645,    0.0204763,   0.00878452,   0.00422461,   0.00210050,   0.00110088],$
    [0.0801816,   0.00916309,   0.00393344,   0.00191854,  0.000959409,  0.000484630],$
    [0.0534391,   0.00609373,   0.00257950,   0.00127993,  0.000651180,  0.000358986],$
    [0.0369379,   0.00421229,   0.00177550,  0.000887714,  0.000452417,  0.000258144]]
sigmadc=[0.0361890,   0.00412712,   0.00173962,  0.000869778,  0.000443276,  0.000252928]
sigmacd=[2.61089e-07,   0.00166342,  0.000794248,  0.000437837,  0.000349688,  0.000313635]

 
 
 data=donnees

 SizeDataMC=size(data)
 SizeImage=SizeDataMC[1]
 SizeMC=SizeDataMC[2]
 image_wt_mc=fltarr(sizeimage,nbrscale1+1,sizemc)
 scaling_wt_mc=fltarr(sizeimage,nbrscale1+1,sizemc)
 vecMC=fltarr(SizeMC)
 trans=fltarr(sizeimage,nbrscale1,sizemc)

 ;coeff_ApAp=fltarr(SizeImage,SizeMC,nbrscale1,nbrscale2)
 ;coeff_AA=fltarr(SizeImage,SizeMC,nbrscale1,nbrscale2)
 ;coeff_AD=fltarr(SizeImage,SizeMC,nbrscale1,nbrscale2)
 ;coeff_DA=fltarr(SizeImage,SizeMC,nbrscale1,nbrscale2)
 ;coeff_DD=fltarr(SizeImage,SizeMC,nbrscale1,nbrscale2)

 support=fltarr(SizeImage,SizeMC,nbrscale1+1,nbrscale2+1)
 support[*]=0

 
 
 cx=fltarr(SizeImage,SizeMC)
 dx=fltarr(SizeImage,SizeMC)
 ccxz=fltarr(SizeImage,SizeMC,nbrscale2+1)
 ccpxz=fltarr(SizeImage,SizeMC,nbrscale2+1)
 cdxz=fltarr(SizeImage,SizeMC)
 dcxz=fltarr(SizeImage,SizeMC)
 ddxz=fltarr(SizeImage,SizeMC)
 temp=fltarr(SizeImage,SizeMC)
 
 cdxztab=fltarr(SizeImage,SizeMC,nbrscale2)
 dcxztab=fltarr(SizeImage,SizeMC,nbrscale1)
 ddxztab=fltarr(SizeImage,SizeMC,nbrscale1,nbrscale2)
 cdxztabp=fltarr(SizeImage,SizeMC,nbrscale2)
 dcxztabp=fltarr(SizeImage,SizeMC,nbrscale1)
 ddxztabp=fltarr(SizeImage,SizeMC,nbrscale1,nbrscale2)

 for j=0,SizeMC-1 do begin
  image=data[*,j]
  mrs_wttrans,image,image_wt,NbrScale=NbrScale1+1
  image_wt_mc[*,*,j]=image_wt.coef
 end
 
 scaling_wt_mc=image_wt_mc
 
 ;for sx=0,nbrscale1-2 do begin
 for sx=0,nbrscale1-1 do begin
  ;;;calcul des coefficients d'échelle de la 1e transformée
  scaling_wt_mc[*,sx,*]=total(image_wt_mc[*,sx:nbrscale1,*],2)
 end

 temp=data
 ccpxz[*,*,0]=temp
 for sz=1,nbrscale2 do begin 
  transformz_approx,temp,cp,cdxyz,0,sz 
  ccpxz[*,*,sz]=cp
  temp=cp
 end
 
 for j=1,nbrscale1 do begin
  print,'scalex =',j
  cx[*,*]=scaling_wt_mc[*,j,*]
  temp=cx
  ccxz[*,*,0]=temp
  

  for k=1,nbrscale2 do begin
   print,'scalez =',k
   ;;calcul des coefficients approx_approx et approx_detail
   transformz_approx,temp,ccx,cdxz,j,k
   ccxz[*,*,k]=ccx
   temp=ccxz[*,*,k]
   ;;Denoising sur les coeff_AD pour j=nbrscale1
   if (j EQ nbrscale1) then begin
    sigma=sigmacd[k-1]
    gaussHardThresholdDenoise,cdxz,alpha,sigma,support,j,k-1
    cdxztab[*,*,k-1]=cdxz
   end
  end
  
  for k=1,nbrscale2 do begin
   ;;calcul des coefficients detail_approx et detail_detail
   ccpxzprev=fltarr(SizeImage,SizeMC)
   ccpxzact=fltarr(SizeImage,SizeMC)
   ccxzprev=fltarr(SizeImage,SizeMC)
   ccxzact=fltarr(SizeImage,SizeMC)
   ccpxzprev[*,*]=ccpxz[*,*,k-1]
   ccpxzact[*,*]=ccpxz[*,*,k]
   ccxzprev[*,*]=ccxz[*,*,k-1]
   ccxzact[*,*]=ccxz[*,*,k]
   transformz_detail,ccpxzprev,ccpxzact,ccxzprev,ccxzact,dcxz,ddxz,j,k
   ;;Denoising sur les coeff_DA pour k=nbrscale2-1
   if (k EQ nbrscale2) then begin
    sigma=sigmadc[j-1]
    gaussHardThresholdDenoise,dcxz,alpha,sigma,support,j-1,k
    dcxztab[*,*,j-1]=dcxz
    end
   ;;Denoising sur les coeff_DD
   sigma=sigmadd[j-1,k-1]
   gaussHardThresholdDenoise,ddxz,alpha,sigma,support,j-1,k-1
   ddxztab[*,*,j-1,k-1]=ddxz
  end
  
 for k=0,nbrscale2 do ccpxz[*,*,k]=ccxz[*,*,k]
 
 end

 


 ;Reconstruction itérative
 if not keyword_set(niter) then niter=5.
 lambda=0.001
 delta=lambda/(niter-1.)
 X=donnees
 X_wt=fltarr(sizeimage,sizemc,nbrscale1)
 tempd_wt_coef=fltarr(sizeimage,sizemc,nbrscale1+1)
 tempd_wt_scaling=fltarr(sizeimage,sizemc,nbrscale1+1)
 X_mswt=fltarr(sizeimage,sizemc,nbrscale1,nbrscale2)
 solution=donnees
 ;solution[*]=0
 solution[*]=1 ;(richardson-lucy)
 solution_wt_coef=fltarr(sizeimage,sizemc,nbrscale1+1)
 solution_wt_scaling=fltarr(sizeimage,sizemc,nbrscale1+1)
 residual_wt_coef=fltarr(sizeimage,sizemc,nbrscale1+1)
 residual_wt_scaling=fltarr(sizeimage,sizemc,nbrscale1+1)
 
 
 chxys=fltarr(sizeimage,sizemc)
 cgxys=fltarr(sizeimage,sizemc,nbrscale1)
 chzs=fltarr(sizeimage,sizemc)
 cgzs=fltarr(sizeimage,sizemc,nbrscale2)
 
 tempd=fltarr(sizeimage,sizemc)
 chxy=fltarr(sizeimage,sizemc)
 cgxy=fltarr(sizeimage,sizemc)
 chz=fltarr(sizeimage,sizemc)
 cgz=fltarr(sizeimage,sizemc)
 
 
 
 
coeff_AA_tab=fltarr(sizeimage,sizemc)
coeff_AD_tab=fltarr(sizeimage,sizemc,nbrscale2)
coeff_DA_tab=fltarr(sizeimage,sizemc,nbrscale1)
coeff_DD_tab=fltarr(sizeimage,sizemc,nbrscale1,nbrscale2)

cxdata=fltarr(sizeimage,sizemc)
dxdata=fltarr(sizeimage,sizemc)
ccxzdata=fltarr(sizeimage,sizemc,nbrscale2+1)
ccpxzdata=fltarr(sizeimage,sizemc,nbrscale2+1)
cdxzdata=fltarr(sizeimage,sizemc)
dcxzdata=fltarr(sizeimage,sizemc)
ddxzdata=fltarr(sizeimage,sizemc)
tempdata=fltarr(sizeimage,sizemc)

cxsol=fltarr(sizeimage,sizemc)
dxsol=fltarr(sizeimage,sizemc)
ccxzsol=fltarr(sizeimage,sizemc,nbrscale2+1)
ccpxzsol=fltarr(sizeimage,sizemc,nbrscale2+1)
cdxzsol=fltarr(sizeimage,sizemc)
dcxzsol=fltarr(sizeimage,sizemc)
ddxzsol=fltarr(sizeimage,sizemc)
tempsol=fltarr(sizeimage,sizemc)



 print,'reconstruction iterative'
 ;tempd=donnees
 for iter=0,niter-1 do begin
 print,'iter = ',iter
 tempd=donnees
 solact=solution
 

 
 if keyword_set(beam) then begin
  for j=0,sizemc-1 do begin
   solutionj=solution[*,j]
   ;if keyword_set(regularization) then beamj=beam[*,j]/(beam[*,j]^2+0.1) else 
   beamj=beam[*,j]
   mrs_convol,solutionj,beamj,solutionjconv
   solution[*,j]=solutionjconv
  end
 endif
 residual=tempd-solution
 
 
 for j=0,sizemc-1 do begin
  residualj=residual[*,j]
  mrs_wttrans,residualj,residual_wt,nbrscale=nbrscale1+1
  residual_wt_coef[*,j,*]=residual_wt.coef
 end
 
 residual_wt_scaling=residual_wt_coef
 ;for sx=1,nbrscale1-1 do begin
 for sx=1,nbrscale1 do begin
  residual_wt_scaling[*,*,sx-1]=total(residual_wt_coef[*,*,sx-1:NbrScale1],3)
 end
 

  
  ;;Decompose XY
  for sx=1,nbrscale1 do begin
   print,'sx=',sx
   
    
   ;;Process cxy_dz if last scalex
   if (sx EQ nbrscale1) then begin
    ;;Decompose Z of cxy (last scalex)
    for sz=1,nbrscale2 do begin
    print,'sz=',sz
     transformz,chxys,chzs,cgzstab,sx,sz
     cgzs[*,*,sz-1]=cgzstab
     
     supportxz=support[*,*,sx,sz-1]
     ;ind=where(supportxz GT 0,taillesup)
     ind=where(supportxz EQ 0,taillesup)
     cgzst=fltarr(sizeimage,sizemc)
     if (taillesup GT 0) then begin
      cgzst[*,*]=cgzs[*,*,sz-1]
      cgzst[ind]=0
      cgzs[*,*,sz-1]=cgzst[*,*]
     end
     cgzst[*,*]=cgzs[*,*,sz-1]
     if not keyword_set(beam) then softthreshold,cgzst,lambda
     cgzs[*,*,sz-1]=cgzst[*,*]
     chxys=chzs
    end
   chxys=chzs+total(cgzs,3)
  endif
   ;;Decompose Z of dxy
   for sz=1,nbrscale2 do begin
   print,'sz=',sz
    cgxyst=fltarr(sizeimage,sizemc)
    cgxyst[*,*]=cgxys[*,*,sx-1]
    transformz,cgxyst,chzs,cgzst,sx,sz
    cgzs[*,*,sz-1]=cgzst[*,*]
     if (sz EQ nbrscale2) then begin
      supportxz=support[*,*,sx-1,sz]
      ;ind=where(supportxz GT 0,taillesup)
      ind=where(supportxz EQ 0,taillesup)
      if (taillesup GT 0) then begin
       ;chzs[ind]=chz[ind]
       chzs[ind]=0
      end
      if not keyword_set(beam) then softthreshold,chzs,lambda
     end
     supportxz=support[*,*,sx-1,sz-1]
     ind=where(supportxz EQ 0,taillesup)
     cgzst=fltarr(sizeimage,sizemc)
     if (taillesup GT 0) then begin
      cgzst[*,*]=cgzs[*,*,sz-1]
      cgzst[ind]=0
      cgzs[*,*,sz-1]=cgzst[*,*]
     end
     cgzst[*,*]=cgzs[*,*,sz-1]
     if not keyword_set(beam) then softthreshold,cgzst,lambda
     cgzs[*,*,sz-1]=cgzst[*,*]
     ;cgxy=chz
     cgxys[*,*,sx-1]=chzs[*,*]
    ;end

   end
   ;;Reconstruct Z of dxy
   cgxys[*,*,sx-1]=chzs+total(cgzs,3)
   solution=chxys
 end
  ;;Reconstruct X
  solution=chxys+total(cgxys,3)

  
  ;Richardson-Lucy
  if keyword_set(beam) then begin
   for j=0,sizemc-1 do begin
    beamj=beam[*,j]
    solactj=solact[*,j]
    solutionj=solution[*,j]
    mrs_convol,solactj,beamj,solactjconv
    denom=solactjconv
    ind=where(denom LT 0.1, taillesup)
    if (taillesup GT 0) then denom[ind]=0.1
    solutionj=(solutionj+solactjconv)/denom
    if keyword_set(regularization) then beamj=beamj/sqrt(beamj^2+0.01)
    mrs_convol,solutionj,beamj,solutionjconv
    solution[*,j]=solutionjconv
   end
   solution=solact*solution
   solution=partie_positive(solution)
  endif else begin
   solution=solact+solution
   solution=partie_positive(solution)
  endelse
  lambda=lambda-delta
  
  
 end


end
