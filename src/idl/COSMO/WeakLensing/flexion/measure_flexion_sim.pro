pro measure_flexion_sim,image,xs,ys,pix_scale,tag=tag,factor=factor,bg=bg,uw=uw

compile_opt idl2
close,101

;; this code takes in a lensed image (from the raytracing suite) and a
;; catalogue of brightness peak positions (from identify_objects.pro),
;; identifies and makes postage stamps of the galaxies, checks for
;; blending and edge effects, deblends images, and measures the
;; flexion on each background galaxy. This is then output to a
;; file. You'll probably need to alter the default file tags
;; below. Inputs that must be specified: image, xs, ys,
;; pix_scale. It is assumed that the input pixel scale is RADIANS PER
;; PIXEL. 


if not keyword_set(bg) then file='lensed' else file='bg'
if not keyword_set(factor) then factor=1.5

nx=(size(image))[1]
ny=(size(image))[2]

openw,101,file+tag+'_flexion.dat'
printf,101,format='(f11.6)',pix_scale*206265.
fmt='(i6,1x,10(f17.6,1x))'
N=(size(xs))[1]

fmax=max(image)
flim=0.02*fmax
for i=0,N-1 do begin
    x_loc=round(xs[i])
    y_loc=round(ys[i])

    loop_centre: xminc=x_loc-10>0
    yminc=y_loc-10>0
    xmaxc=x_loc+10<(nx-1)
    ymaxc=y_loc+10<(nx-1)
    
    cutout=image[xminc:xmaxc,yminc:ymaxc]
    x_loc_new=round(median(where(cutout eq max(cutout)) mod (size(cutout))[1]+$
      xminc))
    y_loc_new=round(median(where(cutout eq max(cutout))/(size(cutout))[1]+$
      yminc))
    
    if x_loc_new ne x_loc or y_loc_new ne y_loc then begin
        x_loc=x_loc_new
        y_loc=y_loc_new
        goto,loop_centre
    endif

    object=search2d(image,x_loc,y_loc,flim,10000)

    xmin=min(object mod nx)
    xmax=max(object mod nx)
    ymin=min(object/nx)
    ymax=max(object/nx)

    if xmin eq 0 or xmax eq nx-1 or ymin eq 0 or ymax eq nx-1 then begin
        message,'object is on the edge',/info
        goto,skip_object
    endif

    range=round(max([xmax-xmin,ymax-ymin])*2.0)
    if range mod 2 eq 1 then range=range+1

    
    if x_loc-range/2 lt 0 then range=2*x_loc
    if x_loc+range/2 gt nx-1 then range=2*(nx-1-x_loc)
    if y_loc-range/2 lt 0 then range=2*y_loc
    if y_loc+range/2 gt nx-1 then range=2*(nx-1-x_loc)
    
    if range gt nx/10 then goto,skip_object

    xmin=x_loc-range/2>0
    xmax=x_loc+range/2<(nx-1)
    ymin=y_loc-range/2>0
    ymax=y_loc+range/2<(nx-1)

    pstamp=image[xmin:xmax,ymin:ymax]

    if (size(pstamp))[1] gt nx/5. or (size(pstamp))[2] gt nx/5. then $
      goto,skip_object

    idx=where(xs gt xmin and xs lt xmax and ys gt ymin and ys lt ymax,count)
    seg=pstamp*0.
    object=search2d(pstamp,x_loc-xmin,y_loc-ymin,flim,10000)
    seg[object]=1
    
    if count gt 1 then begin
        for iobj=0,count-1 do begin
            if idx[iobj] eq i then goto,unblend
            xobj=round(xs[idx[iobj]]-xmin)
            yobj=round(ys[idx[iobj]]-ymin)
            niter=0
            retry: if niter gt 20 then goto,unblend
            if pstamp[xobj,yobj] lt flim then begin
                cutout=pstamp[xobj-5>0:xobj+5<(range-1),$
                              yobj-5>0:yobj+5<(range-1)]
                xobjnew=round(median(where(cutout eq max(cutout)) mod $
                                     (size(cutout))[1]))
                yobjnew=round(median(where(cutout eq max(cutout))/$
                                     (size(cutout))[1]))

                if xobjnew ne xobj or yobjnew ne yobj then begin
                    xobj=xobjnew
                    yobj=yobjnew
                    niter=niter+1
                    goto,retry
                endif
            endif

            if pstamp[xobj,yobj] lt flim then goto,unblend
            newobj=search2d(pstamp,xobj,yobj,flim,10000)
            overlap=where(seg[newobj] gt 0,blendnum)
            if blendnum gt 0 then begin
                message,'object blended',/info
                goto,skip_object
            endif else pstamp[newobj]=min(pstamp)
            unblend:continue
        endfor
    endif

    other=where(pstamp gt flim and seg lt 1,count)
    if count gt 0 then pstamp[other]=min(pstamp)

    e=quickmoments(pstamp,a=a,xc=centroid)
                                ;print,e,a

    if keyword_set(er) then begin
        flexion_er,pstamp,flex,dum1,dum2
        f1=real_part(flex)/pix_scale
        f2=imaginary(flex)/pix_scale
        xc=centroid[0]
        yc=centroid[1]
    
    endif else if keyword_set(uw) then begin
        flexion_moments_uw,pstamp,e,F,G,xc,yc
        f1=F[0]/pix_scale
        f2=F[1]/pix_scale
    endif else begin
        
        
        rhalf=halflight(pstamp)
        weightfunc,pstamp,factor,w,sig2,rhalf=rhalf
        flexion_moments,pstamp,e,F,G,w,sig2,xc,yc


        F=F/pix_scale           ; in rad^(-1)
        G=G/pix_scale           ; in rad^(-1)
        f1=F[0]
        f2=F[1]
    endelse

    xc=(xc+xmin);*pix_scale
    yc=(yc+ymin);*pix_scale

    printf,101,format=fmt,i,xs[i],$
           ys[i],xc,yc,a,e,f1,f2

    ;print,'Schneider-Er says: ',real_part(flex),imaginary(flex)
    ;print,'Goldberg/Leonard says: ',F

    skip_object: continue
endfor
;read,idum

finish:close,101

return
end
