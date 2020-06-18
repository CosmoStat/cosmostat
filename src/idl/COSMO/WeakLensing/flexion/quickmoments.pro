function quickmoments,img,xc=centroid,trace=trace,size=size,a=a,b=b

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

centroid=[xc,yc]

;print,'Flux=',flux
;print,'Center = (',xc,',',yc,')'

Q11=total((i0-xc)^2*img)/flux
Q12=total((i0-xc)*(j0-yc)*img)/flux
Q22=total((j0-yc)^2*img)/flux

trace=Q11+Q22

a=sqrt((Q11+Q22)/2. + sqrt(((Q11-Q22)/2)^2+Q12^2))
b=sqrt((Q11+Q22)/2. - sqrt(((Q11-Q22)/2)^2+Q12^2))

e1=(Q11-Q22)/(Q11+Q22+2.*sqrt(Q11*Q22-Q12^2))
e2=2*Q12/(Q11+Q22+2.*sqrt(Q11*Q22-Q12^2))
size=sqrt(Q11+Q22)

;print,'Size=',size
;print,'Ellipticity=',e1,e2

ell=[e1,e2]
return,ell
end
