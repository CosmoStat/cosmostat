pro weightfunc,img,factor,w,sig2,rhalf=rhalf,mask=mask

;; This program defines a circular gaussian weighting function for a galaxy
;; postage stamp (img). The variance is taken to be some user-specified factor
;; multiplied by the the square of the half-light radius of the image
;; (as determined by halflight.pro). Goldberg & Leonard (2007) found the
;; optimal factor to be around 1.5. The code is also set up to use the
;; sig^2=factor*(Q11+Q22), the quadrupole term approximately
;; describing the effective size of the image. 


nx=(size(img))[1]
ny=(size(img))[2]

usequad=0
if not keyword_set(rhalf) then usequad=1

i0=lonarr(nx*nx)
j0=lonarr(nx*nx)
for i=long(0),nx*nx-1 do begin
    i0[i]=i mod nx
    j0[i]=i/nx
endfor

; First, start with an unweighted measure
flux=total(img)
xc=total(i0*img)/flux
yc=total(j0*img)/flux
Q11=total((i0-xc)^2*img)/flux
Q22=total((j0-yc)^2*img)/flux

niter=3
w=dblarr(nx,nx)
for i=0,niter-1 do begin
     w=dblarr(nx,nx)
     if (usequad eq 1) then begin
         sig2=factor*(Q11+Q22)
     endif else begin 
         sig2=factor*rhalf^2
     endelse
     for i=0,nx-1 do begin
        for j=0,nx-1 do begin
            w[i,j]=exp(-((i-xc)^2+(j-yc)^2)/(2.*sig2))
        endfor
    endfor
    flux=total(w*img)
    xc=total(i0*w*img)/flux
    yc=total(j0*w*img)/flux
    if (usequad eq 1) then begin
        Q11=total((i0-xc)^2*w*img)/flux
        Q22=total((j0-yc)^2*w*img)/flux
    endif
endfor
w=w/total(w)

end
