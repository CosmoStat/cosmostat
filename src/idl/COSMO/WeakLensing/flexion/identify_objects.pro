pro identify_objects,image,fwhm,xs,ys,flim=flim,seg=seg

;; This script identifies the locations of brightness peaks within an
;; image. The inputs are the 2d image and the fwhm of the galaxies
;; generated. This should be the same as used in the gen_gal programme
;; to generate the galaxy images. The outputs xs and ys specify the x
;; and y positions of the sources. 

if n_elements(flim) eq 0 then flim=0.5

;; first find images

find,image,x,y,flux,sharp,roundness,0.5,fwhm,[-2.0,2.0]

N=(size(x))[1]
testseg=image*0
seg=image*0
nim=(size(image))[1]


ns=0
for i=0,N-1 do begin
;    print,i,ns
    if image[round(x[i]),round(y[i])] gt flim*max(image) then begin
        object=search2d(image,round(x[i]),round(y[i]),flim*max(image),max(image)*100.)
        testseg[object]=1.
        ;; identify any other objects
        if i lt N-1 then begin
            xtest=(x[i+1:N-1])
            ytest=(y[i+1:N-1])  ;
            
            idx=where(abs(xtest-x[i]) lt 1.0 and abs(ytest-y[i]) lt 1.0,count)
            if count gt 0 then goto,reloop
            
            xmin=min(object mod nim)
            ymin=min(object/nim)
            xmax=max(object mod nim)
            ymax=max(object/nim)
            idx=where(xtest gt xmin and xtest lt xmax and ytest gt ymin and $
                      ytest lt ymax, count)
            if count gt 0 then goto,reloop
            
        endif
        
        add_source: seg[object]=i+1
        if ns eq 0 then begin
            xs=x[i]
            ys=y[i]
        endif else begin
            xs=[xs,x[i]]
            ys=[ys,y[i]]
        endelse
        ns=ns+1
        reloop: testseg=testseg*0.
    endif
endfor

;; now identify any missed objects in image
;ns=0

reiter:idx=where(image gt flim*max(image) and seg lt 0.1,count)

if count gt 0 then begin
    loc_peak=where(image[idx] eq max(image[idx]))
    x_peak=idx[loc_peak] mod nim
    y_peak=idx[loc_peak]/nim

    object=search2d(image,x_peak,y_peak,flim*max(image),max(image)*100.)
    seg[object]=n_elements(xs)+1
    if ns eq 0 then begin
        xs=x_peak
        ys=y_peak
    endif else begin
        xs=[xs,x_peak]
        ys=[ys,y_peak]
    endelse
    ns=ns+1
;    print,ns
    goto,reiter
endif

return
end
