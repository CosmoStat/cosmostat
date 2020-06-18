;+
; NAME: 
;       mrs_find_obj
;
; PURPOSE: 
;       Detect the sources in  an image by using the Sextractor program 
;       This routine is calling the binary  {sex}. 
;
; CALLING:
;       PRO mrs_find_obj, imag, tabob, tabobj=tabobj, NbrObj=NbrObj, level=level, surf=surf, posneg=posneg, $
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
;                Face            int              Face number
;                Theta           float            Theta coordinate
;                Phi             float            Phi Coordinate
;                SGN             INT              Sig +1 if the structure is positive, -1 otherwise
;                MAXVAL          FLOAT            Maximum value of the object
;                FLUX            FLOAT            Pixel values integration inside the object
;                FWHM            FLOAT            Full width at half maximum
;                AXIS_A          FLOAT            Axis A of the ellipse
;                AXIS_B          FLOAT            Axis B of the ellipse
;                ANGLE           FLOAT            Orientation
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


pro mrs_ell, map, TabObj, maxval=maxval, ImaEll=ImaEll, plot=plot, zoom=zoom
nside = gnside(map)
if keyword_set(zoom) then begin
   nsideout= nside*zoom
    imag = mrs_resize( map, nside=nsideout)   
end else begin
  imag=map 
  zoom=1.
end

vs = size(imag) 
npix = vs[1]
nside = npix2nside(npix)
 
if keyword_set(maxval) then begin
  ind = where (imag GT maxval, c)
  if c gt 0 then imag[ind] = maxval
  ind = where (imag LT -maxval, c)
  if c gt 0 then imag[ind] = -maxval
end
DPatch = TabObj[0].DetectInPatch 
if DPatch EQ 0 then begin
  Face = H2F(Imag)
  FaceEllipse = Face
  FaceEllipse[*]=0
  OneFace = fltarr(nside,nside)
  TabX = TabObj.x
  TabY = TabObj.y
  TabZ = TabObj.Face
  SizeX = nside
  SizeY = nside
end else begin
  Face = mrs_split(Imag)
  FaceEllipse = Face.map
  FaceEllipse[*]=0
  OneFace = fltarr(Face.nx, Face.ny)
  TabX = TabObj.Px
  TabY = TabObj.Py
  TabZ = TabObj.PZ
  SizeX = Face.nx
  SizeY = Face.ny
end
 
vs = size(TabObj)
NbrObj = vs[1]
; print, NbrObj
Ind=0L

MakeImaEll=0
if MakeImaEll EQ 1 then begin

for Ind=0,NbrObj-1 do begin
     OneFace[*]=0
     Obj = TabObj[Ind]
     ;pixf2ang,nside,long(Obj.x+0.5),long(Obj.y+0.5),Obj.face,theta,phi
     ;Obj.theta = theta
     ;Obj.phi = phi
     ;TabObj[Ind] = Obj
    ;if (Obj.axis_a LT 2.5*Obj.axis_b) then begin
       sigx =  Obj.axis_a  / sqrt(2. * alog(2.)) * 3. / 2.
       sigy =  Obj.axis_b  / sqrt(2. * alog(2.)) * 3. / 2.
       ; plot an ellipse at 2sigma
       ;ang2pix,nside,Obj.theta, obj.phi, x, y, f1
       ;if (abs(Obj.x -x) GT 1 or  abs(Obj.y -y) GT 1  or  abs(Obj.face  - f1) GT 1 ) then begin
       ;  print, "Error: ", x,y, f1  , " == ", Obj.x, Obj.y, Obj.face
       ;end
      
       genellipse, TabX[ind]*zoom, TabY[ind]*zoom, Obj.angle, sigx*zoom, sigy*zoom, XEl, YEl
       indel = where( Xel lt SizeX and Xel ge 0 and yel lt SizeY and yel ge 0 and Xel NE TabX[ind] and  Xel NE TabY[ind], c)
       if c GT 0 then  OneFace[XEl[indel], YEl[indel]] = 1
       FaceEllipse[*,*,TabZ[ind]] = FaceEllipse[*,*,TabZ[ind]] + OneFace
    ;end
end

  if DPatch EQ 0 then ImaEll = f2h(FaceEllipse) $
  else begin
     Face.map = FaceEllipse
     ImaEll = mrs_invsplit(Face)
  end
  M = max(imag)
  if keyword_set(maxval) then M = maxval
  ind = where(ImaEll ne 0, c)
  if c GT 0 then  imag[ind] = M
  ; tvs, ImaEll
  if keyword_set(plot) then tvs, imag
  ImaEll = imag
end

DONE:


end


;==================================================

PRO mrs_find_obj, imag,  TabObj, NbrObj=NbrObj, level=level, surf=surf, posneg=posneg, $
  MaxRadiusDetection=MaxRadiusDetection, minfwhm=minfwhm, minmax=minmax, MaxSize=MaxSize, $
  MaxRatioAxis=MaxRatioAxis, imaEll=imaEll, plot=plot, optmr=optmr, maprms=mapmrs, i2=i2, Patch=Patch, frac=frac, SizePatchDegrees=SizePatchDegrees, $
  PTrans=PTrans, mask=mask, param=param, optwtfil=optwtfil, border=border, debugtrue=debugtrue, FaceDebug=FaceDebug, Verb=Verb

if N_PARAMS() LT 2 then begin 
        print, 'CALL SEQUENCE: find_obj, imag,  TabObj, NbrObj=NbrObj, level=level, surf=surf, posneg=posneg, MaxRadiusDetection=MaxRadiusDetection, minfwhm=minfwhm, minmax=minmax, MaxSize=MaxSize, plot=plot'
        goto, DONE
        end 

if not keyword_set(Opt) then Opt = ' '  
if not keyword_set(level) then level=0.2  ; detection level
if not keyword_set(surf) then surf=3     ; minimum number of pixels for a detection

if keyword_set(debugtrue) then begin
disp, win=0
disp, win=1
disp, win=2
disp, win=3
end

NbrTot = 0
vs = size(imag)
npix = vs[1]
nside = npix2nside(npix)
FaceMask = 0
fmask=0
fi2 = 0
if not keyword_set(Patch) then begin
  Face = H2F(Imag)
  if keyword_set(i2) then Facei2 = H2F(i2)
  if keyword_set(maprms) then FaceRMS = H2F(mapmrs)
  if keyword_set(mask) then FaceMask = H2F(mask)
  if keyword_set(debugtrue) then Patchdebug = H2F(debugtrue)
end else begin
 PTrans = mrs_split(Imag, SizePatchDegrees=SizePatchDegrees, frac=frac)
 Face = PTrans.map
 if keyword_set(i2) then begin Facei2 = mrs_split(i2, SizePatchDegrees=SizePatchDegrees, frac=frac) & Facei2 = Facei2.map & end
 if keyword_set(maprms) then begin FaceRMS = mrs_split(mapmrs, SizePatchDegrees=SizePatchDegrees, frac=frac) & FaceRMS = FaceRMS.map & end
 if keyword_set(mask) then begin  FaceMask =  mrs_split(mask, SizePatchDegrees=SizePatchDegrees, frac=frac) & FaceMask = FaceMask.map & end
 if keyword_set(debugtrue) and not keyword_set(FaceDebug) then begin  FaceDebug =  mrs_split(debugtrue, SizePatchDegrees=SizePatchDegrees, frac=frac) & FaceDebug = FaceDebug.map & end
end

if keyword_set(param) then spawn, "\cp -f $ISAP/param/sexparam/*.* ."

vs= size(Face)
Nframe = vs[3]
if keyword_set(Patch) then print, "Number of patches = ", Nframe
TabNbrObjPerFace = fltarr(Nframe)

for i=0,Nframe-1 do begin
   if keyword_set(Verb) then begin
      if not keyword_set(Patch) then print, "Detection on healpix face ", i+1 $
      else print, "Detection on pathch ", i+1
   end
   f = Face[*,*,i]
   if keyword_set(mask) then fmask = FaceMask[*,*,i]
   if keyword_set(i2) then fi2 = Facei2[*,*,i]
   if keyword_set(maprms) then writefits, "weight.fits", FaceRMS[*,*,i]
   if not keyword_set(optmr) then begin
      if keyword_set(optwtfil) then begin
          fi2 = f
	  f[*] = 0
          mr_filter, fi2, f, opt=optwtfil
       end
      find_obj, f,  TabObj, NbrObj=NbrObj, level=level, surf=surf, posneg=posneg, $
           MaxRadiusDetection=MaxRadiusDetection, minfwhm=minfwhm, minmax=minmax, MaxSize=MaxSize, MaxRatioAxis=MaxRatioAxis, i2=fi2, mask=fmask, border=border, Verb=Verb
      if keyword_set(debugtrue) and keyword_set(plot) then begin
        wset, 0
        clear
	wset, 1
        clear
	wset, 2
        clear
	wset, 3
        clear
	wset, 0
	vs = size(f)
	z = 512. / float(vs[1])
        if keyword_set(debugtrue) then load, FaceDebug[*,*,i]
	if keyword_set(debugtrue) then if NbrObj GT 0 then plot_ell, TabObj, zoom=z  
	wset, 1
	load, PTrans.map[*,*,i]
	if NbrObj GT 0 then plot_ell, TabObj, zoom=z  
	wset, 2
	load, f
	if NbrObj GT 0 then plot_ell, TabObj, zoom=z  
	wset, 3
	if keyword_set(mask) then load, fmask
      end
      
      if keyword_set(optwtfil) then PTrans.map[*,*,i] = f
  end else begin 
     mr_find_obj, f, TabObj, NbrObj=NbrObj, level=level, surf=surf, posneg=posneg, $
     MaxRadiusDetection=MaxRadiusDetection, minfwhm=minfwhm, minmax=minmax, MaxSize=MaxSize, MaxRatioAxis=MaxRatioAxis, $
     opt=optmr, mask=fmask, border=border
     Face[*,*,i] = f
  end
     
  TabNbrObjPerFace[i] = NbrObj
  my_command = 'Face'+strcompress(string(i+1), /remove_all) + ' = TabObj '
  ACK = EXECUTE( my_command) 
  NbrTot = NbrTot + NbrObj
  if not keyword_set(Patch) then print, "Detection on healpix FACE ", i+1, " NbrObj = ", NbrObj $
  else print, "Detection on on patch ", i+1, " NbrObj = ", NbrObj
  ; print, TabObj.ELLIPTICITY
end
NbrObj = NbrTot

MyStruc  = { Number : long(0), $
             DetectInPatch: 0, $
             X: 0., $
	     Y: 0., $
	     nside: nside, $
	     Face: 0, $
	     Px: 0., $
	     Py: 0., $
	     Pz: 0., $
	     Theta: 0., $
	     Phi: 0., $
	     Sgn: 1, $
	     MaxVal: 1., $
	     Flux: 1., $
	     ErrorFlux: 0., $
             Fwhm: 0., $
	     axis_a: 0., $
	     axis_b: 0., $
	     Size: 0., $
	     ErrSize: 0., $
	     angle: 0., $
	     Area: 0., $
	     ELLIPTICITY: 0.}
if NbrObj GT 0 then begin
	     
TabObj = replicate(MyStruc, NbrTot )
Ind=0L
print, " TOTAL Number of Objects = ", NbrObj
for i=0,Nframe-1 do begin
  my_command = 'TabFace = Face'+strcompress(string(i+1), /remove_all)  
  ACK = EXECUTE( my_command) 
  for j=0,TabNbrObjPerFace[i]-1 do begin
     Obj = MyStruc
     Obj.nside = nside
     Obj.number = Ind+1
     ObjFace = TabFace[j]
     ; hs, ObjFace
     if not keyword_set(Patch) then begin
	     Obj.x = ObjFace.x
         Obj.y = ObjFace.y
         pixf2ang,nside,long(Obj.x+0.5),long(Obj.y+0.5),i,theta,phi
	    Obj.Face = i
	    Obj.DetectInPatch = 0
     end else begin
       Obj.DetectInPatch = 1
       if j EQ 0 then map_hd_Frame = 0
       LB = mrs_xyz2lb(PTrans, ObjFace.x, ObjFace.y, i) ; , map_hd_Frame=map_hd_Frame)
       ; print,  Obj.x, Obj.y, ' ==> ', lb[0], lb[1]
       lb2ang, lb[0,0], lb[0,1], theta, phi
       ang2pix,nside,theta,phi,x,y,face
       Obj.Face = face
       Obj.x = x
       Obj.y = y
       Obj.Px = ObjFace.x
       Obj.Py = ObjFace.y
       Obj.Pz = i
     end
     Obj.theta = Theta
     Obj.Phi = phi
     Obj.sgn = ObjFace.sgn 
     Obj.MaxVal = ObjFace.MaxVal
     Obj.Flux = ObjFace.Flux
     Obj.ErrorFlux = ObjFace.ErrorFlux
     Obj.Fwhm = ObjFace.Fwhm
     Obj.axis_a = ObjFace.axis_a
     Obj.axis_b = ObjFace.axis_b
     Obj.angle = ObjFace.angle
     Obj.Area = ObjFace.Area
     Obj.ELLIPTICITY= ObjFace.ELLIPTICITY
     TabObj[Ind] = Obj
     ;genellipse, Obj.x, Obj.y, Obj.angle, Obj.axis_a, Obj.axis_b, XEl, YEl
     ; indel = where( Xel lt nside and Xel ge 0 and yel lt nside and yel ge 0, c)
     ; if c GT 0 then  OneFace[XEl[indel], YEl[indel]] = 1
     ; print, "Obj ", Obj.number, " theta = ", Obj.theta, " Phi = ", Obj.Phi, " Fwhm = ", Obj.Fwhm
     Ind = Ind+1
  end
; FaceEllipse[*,*,i] = OneFace
end

; mrs_ell, imag, TabObj, ImaEll=ImaEll
; M = max(imag)
; if keyword_set(plot) then begin
; tvs, imag+ImaEll*M*2
; tvs, ImaEll
; end
end

if keyword_set(param) then begin
  delete, "fpar.sex"
  delete, "default.sex"
  delete, "default.nnw"
  delete, "catparams.param"
  delete, "bgd.fits"
end
; genellipse, 100,120, 30, 20, 10, x, y

DONE:

end

;==================================================

pro obj_to_ascii, filename, TabObj, NbrObj
openw, Unit, filename, /get_lun
printf,  Unit, "ObjNumber   Theta Phi MaxVal  Flux  ErrorFlux FwhmA  FwhmB  "


 for i=0, fix(NbrObj)-1 do begin
    printf, Unit, i+1,  TabObj[i].Theta,  TabObj[i].Phi, TabObj[i].MaxVal, TabObj[i].Flux, TabObj[i].ErrorFlux, TabObj[i].axis_a, TabObj[i].axis_b
 end
 
close, Unit
free_lun, Unit
end

;==================================================

pro obj_to_lb_ascii, filename, TabObj, NbrObj, ring=ring
openw, Unit, filename, /get_lun

if not keyword_set(ring) then printf,  Unit, "l  b  Flux  ErrorFlux"

vs = size(TabObj)
NbrObj = vs[1]
ang2lb, TabObj.theta, TabObj.phi, l ,b
for i=0, long(NbrObj)-1 do begin
    if not keyword_set(ring) then printf, Unit, l[i],  b[i],  TabObj[i].Flux,  TabObj[i].ErrorFlux $
    else begin
      pixf2pix, TabObj[i].nside,TabObj[i].x, TabObj[i].y,TabObj[i].face, ipix   ; ipix nested
      nest2ring, TabObj[i].nside, ipix, pix
      printf, Unit, pix[0], l[i],  b[i],  TabObj[i].Flux,  TabObj[i].ErrorFlux
    end
 end
 
close, Unit
free_lun, Unit
end

;==================================================

pro obj_mask,  TabObj, Map, fwhm=fwhm, nbeam=nbeam, nside=nside
;  Fwhm in arc min
; 
if not keyword_set(nside) then nside = TabObj[0].nside
if not keyword_set(fwhm) then fwhm = 5.
if not keyword_set(nbeam) then nbeam = 3.
beam=fwhm/60./180.*!PI/sqrt(2.*alog(2.))/2.
; print, beam

Theta = TabObj.theta
Phi = TabObj.Phi
ang2pix_nest, nside, Theta, Phi, poslist

Map=fltarr(12.*nside*nside)
Map[*]=1.
for i=0L,n_elements(poslist)-1 do begin
 	pix2vec_nest, nside, poslist[i], vector
        query_disc,nside,vector,nbeam*beam,listpix, /nest
	Map[listpix]=0.
end
; mollview, /on, map
;vs = size(Map)
;nx = vs[1]
;ir = lindgen(nx)
;RING2NEST, Nside, ir, in
;Mapn = Map 
;Mapn[in] = Map
; mrs_tv, mapn
end

;==================================================

pro testm
; ps_mask,30.,1000.,nside=128.,n_sigma_t=5
nside=128L
f = fltarr(nside,nside,64)
for i=0,11 do f[*,*,i] = mygauss(nside,nside, 1)
n = randomn(seed, nside^2*12)*0.01
h = f2h(f)
info, h
hn = h + n
mrs_find_obj, hn*10,  TabObj, NbrObj=NbrObj, /plot, level=1
obj_mask,  TabObj, Map, fwhm=33
tvs, map
end


;==================================================



