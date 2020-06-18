pro flexion_moments_uw,img_in,gam,F,G,xc,yc,a=a2;w,sig2,xc,yc

;; Following Okura et al. (2006) with corrections

;; This code computes the unweighted moments of the galaxy image
;; (img_in) and outputs gam, F, G, xc, yc. If using the ray-tracing
;; code, I've found that even in the absence of noise it's better to
;; use the weighted moments for the calculation, but this code is a
;; useful quick check for simple low-flexion images. 

img=img_in;*w

;; First, compute the ordinary moments.
nx=(size(img))[1]
ny=(size(img))[2]


i0=lonarr(nx*ny)
j0=lonarr(nx*ny)
for i=long(0),nx*ny-1 do begin
    i0[i]=i mod nx
    j0[i]=i/nx
endfor
flux=total(img)
xc=total(i0*img)/flux
yc=total(j0*img)/flux


Q11=total((i0-xc)^2*img)/flux
Q12=total((i0-xc)*(j0-yc)*img)/flux
Q22=total((j0-yc)^2*img)/flux
a2=sqrt(Q11+Q22)

Q111=total((i0-xc)^3*img)/flux
Q112=total((i0-xc)^2*(j0-yc)*img)/flux
Q122=total((i0-xc)*(j0-yc)^2*img)/flux
Q222=total((j0-yc)^3*img)/flux

Q1111=total((i0-xc)^4*img)/flux
Q1112=total((i0-xc)^3*(j0-yc)*img)/flux
Q1122=total((i0-xc)^2*(j0-yc)^2*img)/flux
Q1222=total((i0-xc)*(j0-yc)^3*img)/flux
Q2222=total((j0-yc)^4*img)/flux

;; Go 6th moments!
Q111111=total((i0-xc)^6*img)/flux
Q111112=total((i0-xc)^5*(j0-yc)*img)/flux
Q111122=total((i0-xc)^4*(j0-yc)^2*img)/flux
Q111222=total((i0-xc)^3*(j0-yc)^3*img)/flux
Q112222=total((i0-xc)^2*(j0-yc)^4*img)/flux
Q122222=total((i0-xc)*(j0-yc)^5*img)/flux
Q222222=total((j0-yc)^6*img)/flux


;;;

;Now, compute some components
A=dblarr(4,4)
B=dblarr(4)
F=dblarr(2)
G=dblarr(2)

gam=dblarr(2)
zeta=dblarr(2)
delta=dblarr(2)
eta=dblarr(2)
lambda=dblarr(2)

; Compute gam
denom=Q11+Q22+2*sqrt(Q11*Q22-Q12^2)
gam[0]=(Q11-Q22)/denom
gam[1]=(2*Q12)/denom

xi=Q1111+2*Q1122+Q2222

zeta[0]=(Q111+Q122)/xi
zeta[1]=(Q112+Q222)/xi
delta[0]=(Q111-3*Q122)/xi
delta[1]=(3*Q112-Q222)/xi

eta[0]=(Q1111-Q2222)/xi
eta[1]=2*(Q1112+Q1222)/xi
lambda[0]=(Q1111-6*Q1122+Q2222)/xi
lambda[1]=4*(Q1112-Q1222)/xi

; these indices are zeta and then delta
B[0]=zeta[0];-2*gam[0]*zeta[0]-2*gam[1]*zeta[1]$
  ;-gam[0]*delta[0]-gam[1]*delta[1]
B[1]=zeta[1];-2*gam[1]*zeta[0]+2*gam[0]*zeta[1]$
  ;+gam[1]*delta[0]-gam[0]*delta[1]

B[2]=delta[0];-3*gam[0]*zeta[0]+3*gam[1]*zeta[1]
B[3]=delta[1];-3*gam[0]*zeta[1]-3*gam[1]*zeta[0]

;; The second index is F1,F2,G1,G2
;Ax=B

A[0,0]=0.25*(9.+8*eta[0])
A[0,1]=0.25*(8.*eta[1])
A[0,2]=0.25*(2*eta[0]+lambda[0])
A[0,3]=0.25*(2*eta[1]+lambda[1])

A[1,0]=0.25*(8*eta[1])
A[1,1]=0.25*(-8*eta[0]+9.)
A[1,2]=0.25*(-2*eta[1]+lambda[1])
A[1,3]=0.25*(2*eta[0]-lambda[0])

A[2,0]=0.25*(10*eta[0]+7*lambda[0])
A[2,1]=0.25*(-10*eta[1]+7*lambda[1])
A[2,2]=0.25*(3.)
A[2,3]=0.

A[3,0]=0.25*(10*eta[1]+7*lambda[1])
A[3,1]=0.25*(10*eta[0]-7*lambda[0])
A[3,2]=0.
A[3,3]=0.25*(3.)

;; Put in the corrections due to finite aperture
;A[0,0]=A[0,0]-(3*Q111111+6*Q111122+3*Q112222)/(4.*xi*sig2)
;A[0,1]=A[0,1]-(3*Q111112+6*Q111222+3*Q122222)/(4.*xi*sig2)
;A[0,2]=A[0,2]-(Q111111-2*Q111122-3*Q112222)/(4.*xi*sig2)
;A[0,3]=A[0,3]-(3*Q111112+2*Q111222-Q122222)/(4.*xi*sig2)

;A[1,0]=A[1,0]-(3*Q111112+6*Q111222+3*Q122222)/(4.*xi*sig2)
;A[1,1]=A[1,1]-(3*Q111122+6*Q112222+3*Q222222)/(4.*xi*sig2)
;A[1,2]=A[1,2]-(Q111112-2*Q111222-6*Q122222)/(4.*xi*sig2)
;A[1,3]=A[1,3]-(3*Q111122+2*Q112222-Q222222)/(4.*xi*sig2)

;A[2,0]=A[2,0]-(3*Q111111-6*Q111122-9*Q112222)/(4.*xi*sig2)
;A[2,1]=A[2,1]-(3*Q111112-6*Q111222-9*Q122222)/(4.*xi*sig2)
;A[2,2]=A[2,2]-(Q111111-6*Q111122+9*Q112222)/(4.*xi*sig2)
;A[2,3]=A[2,3]-(3*Q111112-10*Q111222+3*Q122222)/(4.*xi*sig2)

;A[3,0]=A[3,0]-(9*Q111112+6*Q111222-3*Q122222)/(4.*xi*sig2)
;A[3,1]=A[3,1]-(9*Q111122+6*Q112222-3*Q222222)/(4.*xi*sig2)
;A[3,2]=A[3,2]-(3*Q111112-10*Q111222+3*Q122222)/(4.*xi*sig2)
;A[3,3]=A[3,3]-(9*Q111122-6*Q112222+Q222222)/(4.*xi*sig2)

; Now, put in all of the corrections due to the shift in centroid
A[0,0]=A[0,0]-(33*Q11^2+14*Q11*Q22+Q22^2+20*Q12^2)/(4.*xi)
A[0,1]=A[0,1]-(32*Q12*Q22+32*Q11*Q12)/(4.*xi)
A[0,2]=A[0,2]-(3*Q11^2-2*Q11*Q22-Q22^2-4*Q12^2)/(4.*xi)
A[0,3]=A[0,3]-(2*Q11*Q12)/xi

A[1,0]=A[1,0]-(32*Q12*Q22+32*Q11*Q12)/(4*xi)
A[1,1]=A[1,1]-(Q11^2+14*Q11*Q22+20*Q12^2+33*Q22^2)/(4.*xi)
A[1,2]=A[1,2]-(-2*Q12*Q22)/xi
A[1,3]=A[1,3]-(Q11^2+4*Q12^2+Q11*Q22-3*Q22^2)/(4.*xi)

A[2,0]=A[2,0]-3*(11*Q11^2-10*Q11*Q22-Q22^2-20*Q12^2)/(4.*xi)
A[2,1]=A[2,1]-3*(8*Q11*Q12-32*Q12*Q22)/(4.*xi)
A[2,2]=A[2,2]-3*(-2*Q11*Q22+Q11^2+Q22^2+4*Q12^2)/(4.*xi)
A[2,3]=A[2,3]-0.

A[3,0]=A[3,0]-3*(32*Q11*Q12-8*Q12*Q22)/(4.*xi)
A[3,1]=A[3,1]-3*(Q11^2+20*Q12^2+10*Q11*Q22-11*Q22^2)/(4.*xi)
A[3,2]=A[3,2]-0.
A[3,3]=A[3,3]-3*(-2*Q11*Q22+Q11^2+Q22^2+4*Q12^2)/(4.*xi)

;; Put in the new corrections for the limited mask

;A[0,0]=A[0,0]-(-3*Q22*Q1122-9*Q11*Q1111-6*Q12*Q1112-9*Q11*Q1122-6*Q12*Q1222-3*Q22*Q1111)/(4.*xi*sig2)
;A[0,1]=A[0,1]-(-3*Q22*Q1112-9*Q11*Q1222-3*Q22*Q1222-6*Q12*Q1122-9*Q11*Q1112-6*Q12*Q2222)/(4.*xi*sig2)
;A[0,2]=A[0,2]-(3*Q22*Q1122-3*Q11*Q1111-2*Q12*Q1112+9*Q11*Q1122+6*Q12*Q1222-Q22*Q1111)/(4.*xi*sig2)
;A[0,3]=A[0,3]-(-6*Q12*Q1122-9*Q11*Q1112+3*Q11*Q1222-3*Q22*Q1112+Q22*Q1222+2*Q12*Q2222)/(4.*xi*sig2)

;A[1,0]=A[1,0]-(-6*Q12*Q1122-3*Q11*Q1112-9*Q22*Q1112-3*Q11*Q1222-9*Q22*Q1222-6*Q12*Q1111)/(4.*xi*sig2)
;A[1,1]=A[1,1]-(-6*Q12*Q1112-3*Q11*Q2222-6*Q12*Q1222-9*Q22*Q1122-3*Q11*Q1122-9*Q22*Q2222)/(4.*xi*sig2)
;A[1,2]=A[1,2]-(6*Q12*Q1122-Q11*Q1112-3*Q22*Q1112+3*Q11*Q1222+9*Q22*Q1222-2*Q12*Q1111)/(4.*xi*sig2)
;A[1,3]=A[1,3]-(-9*Q22*Q1122-3*Q11*Q1122+Q11*Q2222-6*Q12*Q1112+2*Q12*Q1222+3*Q22*Q2222)/(4.*xi*sig2)

;A[2,0]=A[2,0]+3*(-3*Q22*Q1122+3*Q11*Q1111-6*Q12*Q1112+3*Q11*Q1122-6*Q12*Q1222-3*Q22*Q1111)/(4.*xi*sig2)
;A[2,1]=A[2,1]+3*(-3*Q22*Q1112+3*Q11*Q1222-3*Q22*Q1222-6*Q12*Q1122+3*Q11*Q1112-6*Q12*Q2222)/(4.*xi*sig2)
;A[2,2]=A[2,2]+3*(3*Q22*Q1122+Q11*Q1111-2*Q12*Q1112-3*Q11*Q1122+6*Q12*Q1222-Q22*Q1111)/(4.*xi*sig2)
;A[2,3]=A[2,3]+3*(-6*Q12*Q1122+3*Q11*Q1112-Q11*Q1222-3*Q22*Q1112+Q22*Q1222+2*Q12*Q2222)/(4.*xi*sig2)

;A[3,0]=A[3,0]-3*(-6*Q12*Q1122-3*Q11*Q1112+3*Q22*Q1112-3*Q11*Q1222+3*Q22*Q1222-6*Q12*Q1111)/(4.*xi*sig2)
;A[3,1]=A[3,1]-3*(-6*Q12*Q1112-3*Q11*Q2222-6*Q12*Q1222+3*Q22*Q1122-3*Q11*Q1122+3*Q22*Q2222)/(4.*xi*sig2)
;A[3,2]=A[3,2]-3*(6*Q12*Q1122-Q11*Q1112+Q22*Q1112+3*Q11*Q1222-3*Q22*Q1222-2*Q12*Q1111)/(4.*xi*sig2)
;A[3,3]=A[3,3]-3*(3*Q22*Q1122-3*Q11*Q1122+Q11*Q2222-6*Q12*Q1112+2*Q12*Q1222-Q22*Q2222)/(4.*xi*sig2)

sol=invert(A)#B

F=[sol[0],sol[1]]
G=[sol[2],sol[3]]

end
