;+
; NAME: 
;       find_obj
;
; PURPOSE: 
;       Detect the sources in  an image by using the Sextractor program 
;       This routine is calling the binary  {sex}. 
;
; CALLING:
;       PRO find_obj, imag, tabob, tabobj=tabobj, NbrObj=NbrObj, level=level, surf=surf, posneg=posneg, $
;                            MaxRadiusDetection=MaxRadiusDetection, minfwhm=minfwhm, minmax=minmax, MaxSize=MaxSize, plot=plot
;
; INPUT:
;       Imag: image 
;
; INPUT KEYWORDS:
;        level: float -- Threshold level for the detection. Default value is 0.02
;        surf: float -- Minimum surface an object detection. Default is 3 pixels. 
;        posneg: Scalar -- if set, then negative object are also detected.   
;        MaxRadiusDetection: scalar -- if set, only objects within a distance of  MaxRadiusDetection from the center are detected.  
;        minfwhm: float -- Minimum Full width at half maximum for an object 
;        minmax: float -- All object with a maximum value smaller that minmax are not considered.
;        MaxSize: float -- maximum size for an object. 
;        plot: scalar -- if set, display the image with ellipses overplotted around the detected sources.
;
; OUTPUT KEYWORDS:
;      NbrObj: number of detected objects
;
; OUTPUTS:
;           tabobj: IDL structure which contains the information about 
;              the objects. For each object, we have
;                NUMBER          LONG             Onject number
;                X               FLOAT            X position
;                Y               FLOAT            Y position
;                SGN             INT              Sig +1 if the structure is positive, -1 otherwise
;                MAXVAL          FLOAT            Maximum value of the object
;                FLUX            FLOAT            Pixel values integration inside the object
;                FWHM            FLOAT            Full width at half maximum
;                AXIS_A          FLOAT            Axis A of the ellipse
;                AXIS_B          FLOAT            Axis B of the ellipse
;                ANGLE           FLOAT            Orientation (in degrees)TabObj[Ind].angle
;                AREA            FLOAT            Surface
;                ELLIPTICITY     INT              Ellipticity
;
; EXTERNAL CALLS:
;           Binary sex of Sectractor software
;
; EXAMPLE:
;           detection in an image with all default options 
;                find_obj, Imag, TabObj
;      
;                find_obj, Imag, TabObj, Level=0.2, Surf=3, MaxRadiusDetection=300, minfwhm=0, minmax=0, MaxSize=30, / posneg, /plot
;
; HISTORY:
;	Written: Jean-Luc Starck 2006
;-

;==================================================

pro call_detect, i, a, level=level, surf=surf, i2=i2, checkima=checkima, NameSexParam=NameSexParam
Nl = (size(i))[2]
Nc = (size(i))[1]

;print, 'CALL DETECT'
;info,i 
NameImag =gettmpfilename()
CatName = 'xx1cat.fits'   ;  
; CatName = gettmpfilename() + '.fits'
writefits, NameImag, i

if not keyword_set(NameSexParam) then NameSexParam = '$ISAP/param/sexparam/fpar.sex'
if not keyword_set(i2) then  com = 'sex ' + NameImag  + '  -c  '  + NameSexParam + ' -CATALOG_NAME   ' +  CatName $
else begin
NameImag2 = gettmpfilename()
writefits, NameImag2, i2
com = 'sex   ' + NameImag + '  ' +  NameImag2  + '  -c  '  + NameSexParam + '  -CATALOG_NAME   ' +  CatName
end

if keyword_set(level) then com = com + ' -DETECT_THRESH ' + STRCOMPRESS(STRING(level), /REMOVE_ALL) 
if keyword_set(surf) then com = com + ' -DETECT_MINAREA ' + STRCOMPRESS(STRING(surf), /REMOVE_ALL) 
; print, com
spawn, com
a = mrdfits(CatName ,1,1, /silent)
; checkima = readfits('bgd.fits')
; help, a

delete, NameImag
delete,  CatName
if  keyword_set(i2) then delete, NameImag2
; delete, "bgd.fits"
end

;==================================================

pro cleanobj, i1, a, cross=cross, levelVisu=LevelVisu, MaxSize=MaxSize, zero=zero, minfwhm=minfwhm, minmax=minmax, MaxRatioAxis=MaxRatioAxis, $
     TabOK=TabOK, plot=plot, MaxRadiusDetection=MaxRadiusDetection, imag=imag, mask=mask, border=border, Verb=Verb

if not keyword_set(LevelVisu) then LevelVisu = 10

vs=size(i1)
nx=vs[1]
ny=vs[2]
; MaxRadiusDetection=nx
; MaxSize=nx

if not keyword_set(MaxRadiusDetection) then  MaxRadiusDetection= MAX([nx,ny])*2
if not keyword_set(MaxSize) then MaxSize=MaxRadiusDetection
if not keyword_set(minfwhm) then  minfwhm=-1.
if not keyword_set(minmax) then  minmax=0
if not keyword_set(MaxRatioAxis) then  MaxRatioAxis=100
if not keyword_set(border) then  border=0


if keyword_set(cross) then begin
   tvaxis, i1 < LevelVisu 
end else if keyword_set(plot) then begin
   tvscl, alog(i1+1)
   tvellipse, MaxRadiusDetection,  MaxRadiusDetection, nx/2, ny/2, color=255, thick=2 
end

vs=size(a)
No = vs[1]
TabOK=lonarr(No)
Debug=0
for d=0, no-1 do begin
   DetectOK=1
   a0=a[d] 
   x=[a0.x_world]
   y=[a0.y_world]
   DistCenter = sqrt((a0.x_world - nx/2.)^2. + (a0.y_world - ny/2.)^2.)
   if  DistCenter GT MaxRadiusDetection then DetectOK = 0
   if keyword_set(Debug) then if  DetectOK EQ 0 then print, "DistCenter"
   if (a0.A_IMAGE GT MaxSize) or (a0.b_IMAGE GT MaxSize)  then DetectOK=0
   if keyword_set(Debug) then if  DetectOK EQ 0 then print, "A_IMAGE"
   if a0.B_IMAGE LE minfwhm  then DetectOK=0
   if keyword_set(Debug) then if  DetectOK EQ 0 then print, "B_IMAGE"
   if a0.FLUX_MAX   LE  minmax  then DetectOK=0 
   if keyword_set(Debug) then if  DetectOK EQ 0 then print, "FLUX_MAX"
   if (a0.A_IMAGE GT MaxRatioAxis*a0.b_IMAGE) then DetectOK=0
   if keyword_set(Debug) then if  DetectOK EQ 0 then print, "MaxRatioAxis"
   if keyword_set(imag) then if (imag[fix(a0.x_world),fix(a0.y_world)] EQ  0) then DetectOK=0
   if keyword_set(Mask) then if mask[fix(a0.x_world),fix(a0.y_world)] EQ  0 then DetectOK=0
   DistBorderX = MIN([a0.x_world, nx-a0.x_world])
   DistBorderY = MIN([a0.y_world, ny-a0.y_world])
   DistBorder = min( [DistBorderX, DistBorderY])
   if DistBorder LT border then DetectOK=0
   ; print, a0.x_world, a0.y_world, DistBorderX, DistBorderY, DistBorder
   if keyword_set(Debug) then if  DetectOK EQ 0 then print, "Center at 0"
   if DetectOK EQ 1 then begin
     if keyword_set(cross) then begin 
         
         oplot, x,y, psym=1, color=255
     end else begin
         ; this is 3sigma
         sigx =  a0.A_IMAGE  / sqrt(2. * alog(2.)) * 3. / 2.
	 sigy =  a0.B_IMAGE  / sqrt(2. * alog(2.)) * 3. / 2.
	 ; plot an ellipse at 2sigma
         if keyword_set(plot) then tvellipse, sigx, sigy, a0.x_world, a0.y_world, a0.THETA_IMAGE, color=255
	 ; oplot, x,y, psym=1, color=255
     end
  end
  TabOK[d]=DetectOK
end
ind = where (TabOK EQ 0,c )
indok = where (TabOK EQ 1, c1 )

if keyword_set(Verb) then print, "Number of cleaned object = ", c, " Nbr of Confirmed detections = ", c1
end

;==================================================
 
pro get_my_info,  ap, am=am, tabokp, tabokm=tabokm, TabObj, NbrObj=NbrObj

if type_code(ap) NE 8 then begin
   Nop = 0 
   tabokp = 0
   Np =0 
   Ng = 0
end else begin 
  vs=size(ap)
  Nop = vs[1]
  Np = total(tabokp)
  Ng = Np
end

if keyword_set(am) then begin
   if type_code(am) NE 8 then begin
   Nom = 0 
   tabokm = 0
   Nm = 0
   end else begin 
   vs=size(am)
   Nom = vs[1]
   Nm = total(tabokm)   
   Ng = Np + Nm 
   end
end
NbrObj=Ng

MyStruc  = { Number : long(0), $
             X: 0d, $
	     Y: 0d, $
	     Sgn: 1, $
	     MaxVal: 1d, $
	     Flux: 1d, $
	     ErrorFlux: 0d, $
         Fwhm: 0d, $
         SNR: 0d, $
	     axis_a: 0d, $
	     axis_b: 0d, $
	     Size: 0d, $
	     ErrSize: 0d, $
	     MF_Flux: 0d, $
	     MF_ErrFlux: 0d, $
	     MF_Size: 0d, $
	     MF_ErrSize: 0d, $
         MF_X: 0d, $
	     MF_Y: 0d, $	     
	     angle: 0., $
	     Area: 0., $
	     ELLIPTICITY: 0.}
	     
if Ng GT 0 then TabObj = replicate(MyStruc, Ng ) else TabObj = 0
Ind=long(0)

for i=long(0),Nop-1 do begin
  if tabokp[i] eq 1 then begin
     TabObj[Ind].number = Ind+1
     TabObj[Ind].x = ap[i].x_world
     TabObj[Ind].y = ap[i].y_world
     TabObj[Ind].sgn = 1 
     TabObj[Ind].MaxVal = ap[i].FLUX_MAX
     TabObj[Ind].Flux = ap[i].FLUX_AUTO
     TabObj[Ind].ErrorFlux = ap[i].FLUXERR_AUTO
     TabObj[Ind].Fwhm = ap[i].FWHM_IMAGE
     TabObj[Ind].axis_a = ap[i].A_IMAGE
     TabObj[Ind].axis_b = ap[i].B_IMAGE
     TabObj[Ind].angle = ap[i].THETA_IMAGE
     TabObj[Ind].Area = ap[i].ISOAREA_IMAGE
     TabObj[Ind].ELLIPTICITY= ap[i].ELLIPTICITY
     Ind = Ind+1
  end
end

if keyword_set(am) then begin
  for i=long(0),Nom-1 do begin
    if tabokm[i] eq 1 then begin
     TabObj[Ind].number = Ind+1
     TabObj[Ind].x = am[i].x_world
     TabObj[Ind].y = am[i].y_world
     TabObj[Ind].sgn = -1 
     TabObj[Ind].MaxVal = am[i].FLUX_MAX
     TabObj[Ind].Flux = am[i].FLUX_AUTO
     TabObj[Ind].ErrorFlux = am[i].FLUXERR_AUTO
     TabObj[Ind].Fwhm = am[i].FWHM_IMAGE
     TabObj[Ind].axis_a = am[i].A_IMAGE
     TabObj[Ind].axis_b = am[i].B_IMAGE
     TabObj[Ind].Size = (am[i].A_IMAGE + am[i].B_IMAGE) / 2.
     TabObj[Ind].angle = am[i].THETA_IMAGE
     TabObj[Ind].Area = am[i].ISOAREA_IMAGE
     TabObj[Ind].ELLIPTICITY= am[i].ELLIPTICITY
     Ind = Ind+1
    end
  end
end
end

;==================================================
; sigma = 0.5 * Fwhm / sqrt (2. * log ((double) 2.))
; Fwhm = 2 * sigma * sqrt (2. * log ((double) 2.))
pro plot_ell, TabObj, color=color, zoom=zoom
if not keyword_set(zoom) then zoom=1.
vs=size(TabObj)
No = vs[1]
; print, zoom
for d=0, no-1 do begin
   a0=TabObj[d] 
   x=a0.x
   y=a0.y
 
   ; this is 4 sigma
   KSigma = 4.
   sigx =  a0.axis_a  / sqrt(2. * alog(2.)) * KSigma / 2.
   sigy =  a0.axis_b  / sqrt(2. * alog(2.)) * KSigma / 2.
  ;  print, d+1, x, y, sigx, sigy, a0.angle
   tvellipse, sigx*zoom, sigy*zoom, x*zoom, y*zoom, a0.angle, color=color
   ; oplot, x, y, psym=1, color=255
end

end


;==================================================
 
pro im_ell, map, TabObj, maxval=maxval, ImaEll=ImaEll, plot=plot
vs=size(map)
nx=vs[1]
ny=vs[2]
imag=map
if keyword_set(maxval) then begin
  ind = where (imag GT maxval, c)
  if c gt 0 then imag[ind] = maxval
  ind = where (imag LT -maxval, c)
  if c gt 0 then imag[ind] = -maxval
end

Face = Imag
FaceEllipse = Face
FaceEllipse[*]=0

OneFace = Face 
vs = size(TabObj)
NbrObj = vs[1]
print, NbrObj
Ind=0L
for Ind=0,NbrObj-1 do begin
     OneFace[*]=0
     Obj = TabObj[Ind]
    ;if (Obj.axis_a LT 2.5*Obj.axis_b) then begin
       sigx =  Obj.axis_a  / sqrt(2. * alog(2.)) * 3. / 2.
       sigy =  Obj.axis_b  / sqrt(2. * alog(2.)) * 3. / 2.
       ; plot an ellipse at 2sigma
       genellipse, Obj.x, Obj.y, Obj.angle, sigx, sigy, XEl, YEl
       indel = where( Xel lt nx and Xel ge 0 and yel lt ny and yel ge 0 and Xel NE Obj.x and  Xel NE Obj.y, c)
       if c GT 0 then  OneFace[XEl[indel], YEl[indel]] = 1
       FaceEllipse[*,*] = FaceEllipse[*,*] + OneFace
    ;end
end
ImaEll = FaceEllipse
M = max(imag)
if keyword_set(maxval) then M = maxval
ind = where(ImaEll ne 0, c)
if c GT 0 then  imag[ind] = M
; tvs, ImaEll
if keyword_set(plot) then load, ImaEll
ImaEll = imag

DONE:


end


;==================================================

pro cat2szcat, TabObj,  catalog, NbrObj=NbrObj
vs = size(TabObj)
NbrObj= vs[1]
 clus = {num:0,x:0,y:0,tc:0.,cy:0.,snr:0d}
catalog  =0
 if NbrObj GT 0 and vs[2] EQ 8 then begin
    catalog = REPLICATE(clus, NbrObj)
    catalog.num = 0
    catalog.x =  TabObj.x
    catalog.y = TabObj.y
    catalog.tc = (TabObj.AXIS_A +  TabObj.AXIS_B)/2d0
    catalog.cy =  TabObj.FLUX
    catalog.snr =  TabObj.snr
 end
end
 
;==================================================

PRO find_obj, imag,  TabObj, NbrObj=NbrObj, level=level, surf=surf, posneg=posneg, MaxRatioAxis=MaxRatioAxis, Verb=Verb, $
  MaxRadiusDetection=MaxRadiusDetection, minfwhm=minfwhm, minmax=minmax, MaxSize=MaxSize, plot=plot, color=color, i2=i2, mask=mask, border=border, param=param, NameSexParam=NameSexParam, noclean=noclean

if N_PARAMS() LT 2 then begin 
        print, 'CALL SEQUENCE: find_obj, imag,  TabObj, NbrObj=NbrObj, level=level, surf=surf, posneg=posneg, MaxRadiusDetection=MaxRadiusDetection, minfwhm=minfwhm, minmax=minmax, MaxSize=MaxSize, plot=plot'
        goto, DONE
        end

if not keyword_set(Opt) then Opt = ' '  
if not keyword_set(level) then begin
  if keyword_set(i2) then level = get_noise(i2) / 10. $
  else level = get_noise(imag) / 10.
   ; level=0.2  ; detection level
  end
if not keyword_set(surf) then surf=1     ; minimum number of pixels for a detection
if not keyword_set(color) then color=255     
if keyword_set(param) then spawn, "\cp -f $MRS/sexparam/*.* ." 

Nl = (size(imag))[2]
Nc = (size(imag))[1]

ip = imag > 0
TabOKP = 0
TabOKM = 0
if keyword_set(i2) then begin
      i2p = i2 > 0
      i2p = i2p
end else i2p=0

if max(ip) GT 0 then BEGIN
  call_detect, ip, ap, level=level, surf=surf, i2=i2p, NameSexParam=NameSexParam
   
    ; now we need to clean the detection list from 
    ; spurious detection, introducing our knowledge  
  if type_code(ap) EQ 8 then $
    if not keyword_set(noclean) then cleanobj, ip, ap, TabOK=TabOKP, MaxRadiusDetection=MaxRadiusDetection, minfwhm=minfwhm, minmax=minmax, MaxSize=MaxSize, MaxRatioAxis=MaxRatioAxis, imag=imag, mask=mask, border=border, Verb=Verb $
    else begin NbrObj = N_ELEMENTS(ap.flux_auto) & TabOKP = intarr(NbrObj) + 1 & end
end else ap = -1

if keyword_set(posneg) then begin
   im = imag < 0
   im = -im
   if keyword_set(i2) then begin
      i2m = i2 < 0
      i2m = -i2m
   end else i2m=0
   if max(im) GT 0 then BEGIN
     call_detect, im, am, level=level, surf=surf, i2=i2m, NameSexParam=NameSexParam
     if type_code(am) EQ 8 then $
      if not keyword_set(noclean) then  cleanobj, im, am, TabOK=TabOKM, MaxRadiusDetection=MaxRadiusDetection, minfwhm=minfwhm, minmax=minmax, MaxSize=MaxSize, MaxRatioAxis=MaxRatioAxis, imag=imag, mask=mask, border=border $
      else begin NbrObj = N_ELEMENTS(am.flux_auto) & TabOKM = intarr(NbrObj) + 1 & end
   end else am = -1
end
 
; Extract the information of interest from the detection list
get_my_info,  ap, am=am, tabokp, tabokm=tabokm, TabObj, NbrObj=NbrObj
 
; print, NbrObj
if keyword_set(plot) and NbrObj GT 0 then begin
   window, xsize=nc, ysize=Nl
   tvscl, alog(abs(imag)+1)
   plot_ell, TabObj
   if keyword_set(MaxRadiusDetection) then tvellipse, MaxRadiusDetection,  MaxRadiusDetection, nc/2, nl/2, color=color, thick=2
end


DONE:

end

;==================================================
; sigma = 0.5 * Fwhm / sqrt (2. * log ((double) 2.))
; Fwhm = 2 * sigma * sqrt (2. * log ((double) 2.))

pro  convert_list, tabobj, ap, NbrObj, ip
MyStruc  = { Number : long(0), $
             x_world: 0., $
	     y_world: 0., $
 	     FLUX_MAX: 1., $
	     FLUX_AUTO: 1., $
	     FLUXERR_AUTO: 1., $
             FWHM_IMAGE: 0., $
	     A_IMAGE: 0., $
	     B_IMAGE: 0., $
	     THETA_IMAGE: 0., $
	     ISOAREA_IMAGE: 0., $
	     Size: 0., $
	     ErrSize: 0., $
	     ELLIPTICITY: 0.}

if NbrObj GT 0 then begin	     
   Ap = replicate(MyStruc, NbrObj ) 
  for i=0, NbrObj-1 do begin
     Ap[i].x_world = TabObj[i].PosX
     Ap[i].y_world = TabObj[i].PosY
     Ap[i].FLUX_MAX = ip[ fix(TabObj[i].PosX+0.5), fix(TabObj[i].PosY)]
     Ap[i].FLUX_AUTO = TabObj[i].Flux
     Ap[i].FLUXERR_AUTO = TabObj[i].ErrorFlux
     Ap[i].FWHM_IMAGE = 2.^(TabObj[i].ScaleObj)
     Ap[i].A_IMAGE = 2 * TabObj[i].SigmaX * sqrt (2. * alog ( 2.))
     Ap[i].B_IMAGE = 2 * TabObj[i].SigmaY * sqrt (2. * alog ( 2.))
     if Ap[i].A_IMAGE LT Ap[i].B_IMAGE then begin
       x = Ap[i].A_IMAGE
       Ap[i].A_IMAGE = Ap[i].B_IMAGE
       Ap[i].B_IMAGE = x
     end
     Ap[i].THETA_IMAGE = TabObj[i].Angle
     Ap[i].ISOAREA_IMAGE = 2.^(TabObj[i].ScaleObj+1)
     Ap[i].ELLIPTICITY=  ABS(Ap[i].A_IMAGE - Ap[i].B_IMAGE) / (Ap[i].A_IMAGE + Ap[i].B_IMAGE)   
     ; hs, Ap[i]
  end
end

end

;==================================================

PRO mr_find_obj, imag,  TabObj, NbrObj=NbrObj, level=level, surf=surf, posneg=posneg, MaxRatioAxis=MaxRatioAxis, $
  MaxRadiusDetection=MaxRadiusDetection, minfwhm=minfwhm, minmax=minmax, MaxSize=MaxSize, mask=mask, plot=plot, color=color, opt=opt, ImaDetect=ImaDetect, border=border

if N_PARAMS() LT 2 then begin 
        print, 'CALL SEQUENCE: mr_find_obj, imag,  TabObj, NbrObj=NbrObj, level=level, surf=surf, posneg=posneg, MaxRadiusDetection=MaxRadiusDetection, minfwhm=minfwhm, minmax=minmax, MaxSize=MaxSize, plot=plot'
        goto, DONE
        end

if not keyword_set(Opt) then Opt = ' '  
if not keyword_set(level) then level=0.2  ; detection level
if not keyword_set(surf) then surf=3     ; minimum number of pixels for a detection
if not keyword_set(color) then color=255     

Nl = (size(imag))[2]
Nc = (size(imag))[1]

ip = imag > 0
TabOKP = 0
TabOKM = 0
ImaDetect = fltarr(Nc, Nl)

if max(ip) GT 0 then BEGIN
  mr_detect, ip, ap, NbrObj=NbrObj, tabobj=tabobj, opt=opt
  ; print, "NbrObj= " , NbrObj
  if NbrObj GT 0 then ImaDetect = ap
  convert_list, tabobj, ap, NbrObj, ip
  ; hs, ap
    ; now we need to clean the detection list from 
    ; spurious detection, introducing our knowledge  
  if type_code(ap) EQ 8 then $
    cleanobj, ip, ap, TabOK=TabOKP, MaxRadiusDetection=MaxRadiusDetection, minfwhm=minfwhm, minmax=minmax, MaxSize=MaxSize, MaxRatioAxis=MaxRatioAxis, mask=mask, border=border
    ; print, TabOKP
end else ap = -1

if keyword_set(posneg) then begin
   im = imag < 0
   im = -im
   if max(im) GT 0 then BEGIN
     mr_detect, im, am, NbrObj=NbrObj, tabobj=tabobj, opt=opt
     convert_list, tabobj, am, NbrObj, ip
     if type_code(am) EQ 8 then $
      cleanobj, im, am, TabOK=TabOKM, MaxRadiusDetection=MaxRadiusDetection, minfwhm=minfwhm, minmax=minmax, MaxSize=MaxSize, MaxRatioAxis=MaxRatioAxis, mask=mask, border=border
    if NbrObj GT 0 then ImaDetect = ImaDetect - am
   end else am = -1
end
 
; Extract the information of interest from the detection list
get_my_info,  ap, am=am, tabokp, tabokm=tabokm, TabObj, NbrObj=NbrObj
 
; print, NbrObj
if keyword_set(plot) and NbrObj GT 0 then begin
   window, xsize=nc, ysize=Nl
   tvscl, alog(abs(imag)+1)
   plot_ell, TabObj
   if keyword_set(MaxRadiusDetection) then tvellipse, MaxRadiusDetection,  MaxRadiusDetection, nc/2, nl/2, color=color, thick=2
end


DONE:

end
