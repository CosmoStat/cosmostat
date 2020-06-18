;+
; NAME:
;        mrs_msvsts_iuwt_param_computing
;
; PURPOSE:
;     For a given number of scales, determines the VST operator at each scale for the MS-VST
;     transform with spherical isotropic Undecimated Wavelet Transform. At scale j, the VST operator is
;     Tj(aj) = b(j) * sgn(aj+c(j)) * sqrt(abs(aj+c(j))) where aj is the jth scale coefficient of the wavelet transform.
;
;
; CALLING:
;
;     mrs_msvsts_iuwt_param_computing,nbr,c,b,h,tau1,tau2,tau3,sigma
;       
;
; INPUTS:
;     nbr -- number of scales for the MS-VST transform
;    
; OUTPUTS:
;     c -- fltarr[nbr] : vector of the c(j) coefficients for each scale j
;     b -- fltarr[nbr] : vector of the b(j) coefficients for each scale j
;     h -- fltarr[size_ondelette,nbr] : h[*,j] is the low pass filter which gives the jth scale from the original image
;     tau1 -- fltarr[nbr] : vector of the 1st order moments of h[*,j] for each scale j
;     tau2 -- fltarr[nbr] : vector of the 2st order moments of h[*,j] for each scale j
;     tau3 -- fltarr[nbr] : vector of the 3rd order moments of h[*,j] for each scale j
;     sigma -- fltarr[nbr-1] : vector of the asymptotic standard deviations of detail coefficients issued from locally homogeneous 
;     parts of a signal for each wavelet scale
;
;
; EXTERNAL CALLS:
;       mrs_wttrans (mrs)
;
; EXAMPLE:
;
;       Compute the VST operator for a MS-VST transform with 6 scales
;               mrs_msvsts_iuwt_param_computing,6,c,b,h,tau1,tau2,tau3,sigma
;         
; HISTORY:
;	Written: Jérémy Schmitt & Jean-Luc Starck, 2009
;	February, 2009 File creation
;--------------------------------------------------------------------------------------------------------


pro mrs_msvsts_iuwt_param_computing,nbr,c,b,h,tau1,tau2,tau3,sigma

n=ulong(512)
im=fltarr(n*n*12)

f=h2f(im)
f[n/2,n/2,4]=1
im=f2h(f)

mrs_wttrans,im,w,nbrscale=nbr

size_ondelette=size(w.coef)
size_ondelette=size_ondelette[1]

h=fltarr(size_ondelette,nbr)

h[*,nbr-1]=w.coef[*,nbr-1]
for j=0,nbr-2 do begin
 h[*,j]=total(w.coef[*,j:nbr-1],2)
end
	
tau1=fltarr(nbr)
tau2=fltarr(nbr)
tau3=fltarr(nbr)

tau1=total(h,1)
tau2=total(h^2,1)
tau3=total(h^3,1)

c=fltarr(nbr)
c=7*tau2/(8*tau1)-tau3/(2*tau2)

b=fltarr(nbr)
b=sgn(tau1)/sqrt(tau1)

sigmacarre=fltarr(nbr-1)
sigma=fltarr(nbr-1)
for j=0,nbr-2 do begin
 sigmacarre[j]=(tau2[j]/(4*tau1[j]^2))+(tau2[j+1]/(4*tau1[j+1]^2))-(total(h[*,j]*h[*,j+1])/(2*tau1[j]*tau1[j+1]))
 sigma[j]=sqrt(sigmacarre[j])
end


save,filename='resultat_calcul_vst.xdr',c,b,h,tau1,tau2,tau3,im,w,sigma

end