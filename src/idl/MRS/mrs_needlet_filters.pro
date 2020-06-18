;+
; NAME:
;        mrs_needlet_filters
;
; PURPOSE:
;	Calculate the needlet filter, which can be used in an undecimated wavelet transform on the sphere (see mrs_wttrans.pro).
;
; CALLING:
;      Filter = mrs_needlet_filters(Lmax, B=B, Nscale=Nscale, plot=plot, tabphi=tabphi, tabpsi=tabpsi, TabFilterH=TabFilterH, TabFilterG= TabFilterG)
;       
; INPUTS:
;     Lmax -- int: maximum value of l in the spherical harmonic decomposition.
;
; INPUT KEYWORDS:
;      B -- float :  needlet parameter. Default is 2.
;      plot -- scalar: if set, the filters are plotted in a window
;
; INPUT/OUTPUT KEYWORDS:
;      Nscale -- int: Number of bands in the decomposition. By default is automatically calculated:
;                         Nscale = log(lmax)/ log(B)
;
; OUTPUT KEYWORDS:
;       Tabphi -- IDL array[0:lmax, Nscale-1]  = Scaling function  at resolution level
;       Tabpsi -- IDL array[0:lmax, Nscale-1]   = Wavelet function  at resolution level
;       TabFilterH -- IDL array[0:lmax, Nscale-1] = filter H allowing to go from a resolution to the next one
;       TabFilterG -- IDL array[0:lmax, Nscale-1] = filter G allowing to compute the wavelet coeff from the previous resolution level
;
; HISTORY:
;	Written: Jeremy Schmitt, Florent Sureau, JL Starck, 2010
;-------------------------------------------------------------------------------------------------------------------------------------


 ;=====================================================

;;Construction de la fonction f(t)
;;f(t) = exp(-1/(1-t^2)) , -1=<t=<1
;;        0                 ailleurs

function nd_f,t
 if t GT -1 and t LT 1 then begin y=exp(-1./(1.-t^2))
 endif else y=0.
 return,y
end
 ;=====================================================

;;Construction de la fonction psi(u)
function nd_psi,u, den=den
if(u lt -0.90) then num=0. else num=qromb('nd_f',-1,u,/double)
 if not keyword_set(den) then den=qromb('nd_f',-1,1,/double)
 y=num/den
 return,y
end
 ;=====================================================

;;Constuction de la fonction phi(t)
function nd_phi,t,B, den=den
 if t GE 0 and t LE 1/B then y=1.
 if t GT 1/B and t LT 1 then y=nd_psi(1.-2.*B*(t-1./B)/(B-1.), den=den)
 if t GE 1 then y=0.
 return,y
end
 ;=====================================================

;;Construction de la fonction filtre b(xi) de support [1/B,B]
;;Ex : y=filter(0,2.) : calcule b(0), b ÔøΩtant la fonction filtre de support [1/2,2]
function nd_filter,xi,B, den=den
 y=sqrt(abs(nd_phi(xi/B,B, den=den)- nd_phi(xi,B, den=den))) ;Possible that phi(xi/B,B)-phi(xi,B)< 0 due to numerical errors
 return,y
end

 ;=====================================================

;;filter_tab : calcule la valeur de la fonction filtre b(xi) sur un tableau de xi
;;Ex : tabx=findgen(101.)/100.*5. (tableau de valeurs entre 0 et 5)
;;     taby=filter_tab(tabx,2.) (tableau des valeurs du filtre pour xi entre 0 et 5)
;;     plot,tabx,taby
function nd_filter_tab,tabx,B, den=den
 sizetab=size(tabx)
 sizetab=sizetab[1]
 taby=tabx
 for j=0,sizetab-1 do begin
  taby[j]=nd_filter(tabx[j],B, den=den)
 end
 return,taby
end

 ;=====================================================

function nd_filter_jtab,tabx, B0, K, den=den
 B = B0^K
 sizetab=size(tabx)
 sizetab=sizetab[1]
 taby=tabx
 minL=ceil(B/B0)
 maxL=floor(B*B0)
 ; print,minL,maxL
 if(minL gt 0) then taby[0:minL]=0. else minL=0; 
 if((sizetab-1) gt maxL) then taby[maxL:sizetab-1]=0. else maxL=sizetab-1;
; print,"from ",minL," to ",maxL," size ",sizetab
 for j=minL,maxL do begin
  taby[j]=nd_filter(tabx[j]/B,B0, den=den)
 end
 if(K eq 0) then taby[0]=1 ;Approximation coefficient 
 return,taby
end

 ;=====================================================

function needlet_get_all_scales, B, Lmax, Nscale=Nscale, plot=plot
tabx = findgen(Lmax+1)
B=double(B)
Lmax=double(Lmax)
Nscale=ceil(alog(lmax)/alog(B))
if(Lmax eq B^(Nscale-1)) then Nscale=Nscale-1 ;b(B)=0
F=dblarr(lmax+1,Nscale+1)
if keyword_set(plot) then window,/free
for k=0, Nscale do begin
	 F[*, Nscale-k] = nd_filter_jtab(tabx,B, k, den=den)
 
 	if keyword_set(plot) then begin
  		  if k EQ 0 then plot, F[*,Nscale-k], xrange=[0,Lmax+10] , yr=[0,1.2] $
  		  else oplot, F[*,Nscale-k]
 	end
end
Nscale= Nscale+1

return,F
end

;=====================================================

function mrs_needlet_filters,  LmaxIn,  B=B, plot=plot, tabphi=tabphi, tabpsi=tabpsi, TabFilterH=TabFilterH, TabFilterG= TabFilterG, Nscale=Nscale

if not keyword_set(B) then B=2.0
  
tabx = findgen(LmaxIn +1)
B=double(B)
Lmax=double(LmaxIn)

if not keyword_set(Nscale) then F = needlet_get_all_scales(B, Lmax, Nscale=Nscale, plot=plot)  $
else BEGIN
	F=dblarr(Lmax+1,Nscale)
	Lastscale=ceil(alog(Lmax)/alog(B))
	if(Lmax eq B^(Lastscale-1)) then Lastscale=Lastscale-1 ;b(B)=0
	if(Nscale gt Lastscale) then  Nscale=Lastscale
  	
	if(Nscale lt 2) then F=dblarr(lmax+1,Nscale)+1. $
	else begin
 			if keyword_set(plot) then window,/free
  	       for k=0,NScale-1 do begin
 	 	  			F[*,k]= nd_filter_jtab(tabx,B, Lastscale-k, den=den)
  		 		    if(k eq NScale-1) then begin
   			   				minL=(B^(Lastscale-k))
   						   F[0:minL,k]=1.
    	           endif
   	              if keyword_set(plot) then begin
      					  if k EQ 0 then plot, F[*,k], xrange=[0,Lmax+10] , yr=[0,1.2] $
     					  else oplot, F[*,k]
   	   			endif
  		 endfor
    endelse
endelse
;help, F
;print, "FIL PSI-PHI", nscale 
    TabFilterH = fltarr( lmax+1, Nscale-1)
    TabFilterG = fltarr(lmax+1, Nscale-1)
    TabPhi = fltarr( lmax+1, Nscale-1)
    TabPsi = fltarr( lmax+1, Nscale-1)
    
    TabPsi = F[*,0:Nscale-2]
    TabPhi[*, Nscale-2] =  F[*,Nscale-1]
    for j=Nscale-3,0,-1 do  TabPhi[*, j] = sqrt(  TabPhi[*, j+1]^2 +  TabPsi[*, j+1]^2 )
    j=0
    TabFilterH[*,0] = TabPhi[*,0]
    TabFilterG[*,0] = TabPsi[*,0]
    for j=1,Nscale-2 do begin
       TabFilterH[*,j] = TabPhi[*,j] / TabPhi[*,j-1]
       TabFilterG[*,j] = TabPsi[*,j] / TabPhi[*,j-1]
    end


return,F
end

;=====================================================

pro needlets,  imin, imout, B=B, nbrscale=nbrscale, rev=rev, Filter=Filter

if not keyword_set(B) then B=2.0
if not keyword_set(nbrscale) then nbrscale=6

if not keyword_set(rev) then begin
 	sizeimage=(size(imin))[1]
 	imout=fltarr(sizeimage,nbrscale)
    mrs_almtrans,imin,imalm, /tab
     if keyword_set(filter) then flt=filter $
     else flt= mrs_needlet_filters(nscale=nbrscale, imalm.lmax, B=B)
  	for n=0,nbrscale-1 do begin
       g = flt[*,n]
       A = imalm
       mrs_alm_convol, A, g
       mrs_almrec, A,  w_coef
       imout[*,n]= w_coef
    endfor
endif else begin
    vs = size(imin)
    Np = vs[1]
    nbrscale = vs[2]
    print, nbrscale
    imout=fltarr(Np)
    mrs_almtrans, imin[*, 0], AlmBand, /tab
    if keyword_set(filter) then flt=filter $
    else flt= mrs_needlet_filters(nscale=nbrscale, AlmBand.lmax, B=B) 
    ASol = AlmBand.Alm
    ASol[*] = 0
    for n=0,nbrscale-1 do begin
     	if n GT 0 then mrs_almtrans, imin[*, n] , AlmBand, /tab
    	g = flt[*,n]
    	mrs_alm_convol, AlmBand, g
    	ASol = ASol + AlmBand.Alm
    endfor
    AlmBand.Alm = ASol
    mrs_almrec, AlmBand,  imout
endelse
Filter = flt
end

;=====================================================

pro test_needlets, im, imrec

im=randomn(seed,128d*128d*12d)
mrs_almtrans,im,imalm
szalm=(size(imalm.alm))[1]
off=(128d*129d)*2d/2d;ensure no coefs after 2*Nside (very bad decomposition/reconstruction otherwise)
print,off,szalm
imalm.alm(off:szalm-1,0)=0.
imalm.alm(off:szalm-1,1)=0.
mrs_almrec, imalm,im2

needlets,im2,coefs,B=2,nbrscale=6,rev=0 ;Compute coefs per scale
needlets,coefs,imrec,B=2,nbrscale=6,rev=1 ;Reconstruct image per scale

imsum=dblarr(128d*128d*12d) ;sum of images per scale
for k=0,5 do begin
  mollview,imrec[*,k],nested=1
  imsum=imsum+imrec[*,k]
endfor

mollview,imsum-im2,nested=1 ;residual

end

;=====================================================


pro cmp_needlet_cubic

i = fltarr(256,256,12)
i[128,128,4] = 100.
h = F2H(i)
mrs_wttrans, h, w, nbrs=-1
mrs_wttrans, h, w1, nbrs=-1, /need
mrs_wttrans, h, w2, nbrs=-1, /meyer

j=3
w.coef[*,j] =  w.coef[*,j] / max(w.coef[*,j])
tvs, w.coef[*,j], title='Spline  Wavelet', png='fig_splinew.png'

w1.coef[*,j+2]  = w1.coef[*,j+2]  / max(w1.coef[*,j+2])
tvs, w1.coef[*,j+2], title='Needlet  Wavelet', png='fig_needletw.png'

w2.coef[*,j] =  w2.coef[*,j] / max(w2.coef[*,j])
tvs, w2.coef[*,j], title='Meyer  Wavelet', png='fig_meyerw.png'

FaceSpline = (h2f(w.coef[*,j]))[*,*,4] 
FaceNeedlet1 = (h2f(w1.coef[*,j+2]))[*,*,4]  
FaceMey = (h2f(w2.coef[*,j]))[*,*,4] 
 
setps, filename='fig_splinew_versus_needletw.ps', /portrait
plot,  (FaceSpline[*,128]), xrange=[0,250], yr=[-0.4, 1.1], xticks=5, title='Spline wavelet (continuous) versus Needlet wavelet (dotted)'
oplot,  (FaceNeedlet1[*,128]), line=2
endps

setps, filename='fig_abssplinew_versus_absneedletw.ps', /portrait
plot,  abs(FaceSpline[*,128]), xrange=[0,250], yr=[0, 0.4], xticks=5, title='ABS(Spline wavelet) (continuous) versus ABS(Needlet wavelet) (dotted)'
oplot,  abs(FaceNeedlet1[*,128]), line=2
endps

setps, filename='fig_meyerw_versus_needletw.ps', /portrait
plot,  (FaceMey[*,128]), xrange=[0,250], yr=[-0.4, 1.1], xticks=5, title='Meyer wavelet (continuous) versus Needlet wavelet (dotted)'
oplot,  (FaceNeedlet1[*,128]), line=2
endps

setps, filename='fig_absmeyerw_versus_absneedletw.ps', /portrait
plot,  abs(FaceMey[*,128]), xrange=[0,250], yr=[0, 0.4], xticks=5, title='ABS(meyer wavelet) (continuous) versus ABS(Needlet wavelet) (dotted)'
oplot,  abs(FaceNeedlet1[*,128]), line=2
endps

end

;=====================================================


