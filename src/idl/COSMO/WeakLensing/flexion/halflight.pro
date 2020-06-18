function halflight,img,x0=x0

;; This function identifies the halflight radius of a galaxy image
;; based on the input postage stamp (img). As an optional output, x0
;; is a vector specifying the location of the image centroid
;; Note: no weighting is used in this calculation.

; First, the as the center the middle of the image.
F=total(img)

nx=(size(img))[1]
ny=(size(img))[2]
xc=nx/2.
yc=nx/2.

i0=lonarr(nx*ny)
j0=lonarr(nx*ny)
for i=long(0),nx*ny-1 do begin
    i0[i]=i mod nx
    j0[i]=i/nx
endfor
flux=total(img)
xc=total(i0*img)/flux
yc=total(j0*img)/flux

dr=0.1
nbins=fix(nx/2./dr)
for iter=0,2 do begin
    for ibin=0,nbins-1 do begin
        r=dr*(ibin+1)
        idx=where((i0-xc)^2+(j0-yc)^2 lt r^2,n)
        if (n gt 0) then begin
            frac=total(img[idx])/F
        endif else begin
            frac=0.
        endelse
        if (frac gt 0.5) then begin
            ibin=nbins
            xc=total(i0[idx]*img[idx])/total(img[idx])
            yc=total(j0[idx]*img[idx])/total(img[idx])
        endif    
    endfor
endfor
x0=[xc,yc]


return,r

end
