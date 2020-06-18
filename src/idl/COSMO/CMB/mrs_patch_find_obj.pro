;+
; NAME: 
;       mrs_patch_find_obj
;
; PURPOSE: 
;       Detect the sources in  an image by using the Sextractor program 
;       This routine is calling the binary  {sex}. 
;
; CALLING:
;       PRO mrs_patch_find_obj, imag, tabob, tabobj=tabobj, NbrObj=NbrObj, level=level, surf=surf, posneg=posneg, $
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

  PRO mrs_patch_find_obj,  PTrans,  TabObj,  NbrObj=NbrObj, level=level, surf=surf, posneg=posneg, $
  MaxRadiusDetection=MaxRadiusDetection, minfwhm=minfwhm, minmax=minmax, MaxSize=MaxSize, $
  MaxRatioAxis=MaxRatioAxis, imaEll=imaEll, Prms=Prms, Pi2=Pi2,  $
  Pmask=Pmask, border=border, debugtrue=debugtrue, PDebug=PDebug, Verb=Verb, ConvFluxFactor=ConvFluxFactor, $
  map_hd=map_hd, szfit=szfit, pixsurf=pixsurf, TabPatch=TabPatch

COMMON C_PLANCK


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
nside =  PTrans.nside
FaceMask = 0
fmask=0
fi2 = 0

if keyword_set(param) then spawn, "\cp -f $MRS/param/sexparam/*.* ."

vs= size(PTrans.map )
Nframe = vs[3]
print, "Number of patches = ", Nframe
TabNbrObjPerFace = fltarr(Nframe)

FromPatch = 0L   
 if keyword_set(szfit) then FromPatch = 114L   ; i.e. All patches with | b | > 20 degrees.
 
for i=FromPatch,Nframe-1 do begin
; for i=300,300  do begin
   if keyword_set(Verb) then  print, "Detection on patch ", i+1
   f =  PTrans.map[*,*,i]
   if keyword_set(Pmask) then fmask = PMask.map[*,*,i]
   if keyword_set(Pi2) then fi2 = Pi2.map[*,*,i]
   if keyword_set(Prms) then writefits, "weight.fits", PRMS.map[*,*,i]
   level=0
   surf=0
   find_obj, f,  TabObj, NbrObj=NbrObj, level=level, surf=surf,  $
           MaxRadiusDetection=MaxRadiusDetection, minfwhm=minfwhm, minmax=minmax, MaxSize=MaxSize, MaxRatioAxis=MaxRatioAxis, i2=fi2, mask=fmask, border=border, Verb=Verb

   ; This part is only SZ cluster detections in PLANCK maps
   if keyword_set(SZFIT) and NbrObj GT 0 then begin
      if keyword_set(ConvFluxFactor) then fi2 = fi2 * ConvFluxFactor
      
      Flux1 = TabObj.Flux
      TabObj.flux = TabObj.flux * ConvFluxFactor *  pixsurf 
      TabObj.Size = TabObj.FWHM * sqrt(pixsurf)
      TabObj.ERRSIZE = TabObj.ERRSIZE  * sqrt(pixsurf)
      TabObj.ErrorFlux = TabObj.ErrorFlux * ConvFluxFactor *  pixsurf
      cat2szcat, TabObj,  detcat
      
      ; If TabPatch is set, the clusters flux are obtained used a multichannel matched filter
      ; while if not, flux are derived on the SZ map using a single channel matched filter
      if keyword_set(TabPatch) then begin
        vs = size(TabPatch)
        FirstFreq=3
        print, " ===> MMF Flux "
        NbrChannels = P_NBRCANNELS - FirstFreq
        TabNu = P_nuGHz[FirstFreq:P_NBRCANNELS-1]
        TabConvert    =  conversionfactor(TabNu, /antenna2thermo) / ( 2.72500 * 1.e3)
        Nx = (size(f))[1]
        Ny = (size(f))[2]
        TabFrame = fltarr(Nx,Ny,NbrChannels)
        for  c=0,NbrChannels-1 do TabFrame[*,*,c] = tabpatch[c].map[*,*,i] * TabConvert[c]
        nc =  planck_mmf_flux(TabFrame, DetCat, map_hd=map_hd)
      end else nc = planck_smf_get_cluster(fi2, detcat, map_hd)
      
 	     ;; MF_ErrFlux: 0d, $
 	     ;; MF_ErrSize: 0d, $
      TabObj.MF_X = nc.x
      TabObj.MF_Y = nc.y
      TabObj.MF_flux = nc.cy
      TabObj.MF_Size = nc.tc
      TabObj.MF_ErrSize = nc.err_tc
      TabObj.MF_ErrFlux = nc.err_cy
      TabObj.SNR = nc.sn
      Flux2 =      TabObj.MF_flux  
      ; print, "Flux = ", Flux1, Flux1 * ConvFluxFactor * pixsurf,  Flux2
   endif  ; SZ case  
   
    if keyword_set(debugtrue) then begin
    vs = size(f)
	z = 512. / float(vs[1])
	wset, 0
     if keyword_set(debugtrue) then load, PDebug.map[*,*,i]
	if keyword_set(debugtrue) then if NbrObj GT 0 then plot_ell, TabObj, zoom=z  
	wset, 1
	load, fi2
    if NbrObj GT 0 then plot_ell, TabObj, zoom=z  
	wset, 2
	load, PTrans.map[*,*,i]
	if NbrObj GT 0 then plot_ell, TabObj, zoom=z  
	wset, 3
	if keyword_set(pmask) then load, fmask
    end
           
  TabNbrObjPerFace[i] = NbrObj
  my_command = 'Face'+strcompress(string(i+1), /remove_all) + ' = TabObj '
  ACK = EXECUTE( my_command) 
  NbrTot = NbrTot + NbrObj
   print, "Detection on on patch ", i+1, " NbrObj = ", NbrObj
  ; print, TabObj.ELLIPTICITY
end
NbrObj = NbrTot

; GOTO,  DONE

MyStruc  = { Number : long(0), $
             DetectInPatch: 0, $
             X: 0., $
	     Y: 0., $
	     nside: nside, $
	     Face: 0, $
	     Px: 0d, $
	     Py: 0d, $
	     Pz: 0d, $
	     Theta: 0d, $
	     Phi: 0d, $
	     Sgn: 1, $
	     MaxVal: 1d, $
	     Flux: 0d, $
	     ErrorFlux: 0d, $
             Fwhm: 0d, $
             SNR: 0d, $
             MF_ErrFlux: 0d, $
 	     MF_ErrSize: 0d, $
             MF_flux: 0d, $
             MF_Size: 0d, $
	     MF_X: 0d, $
	     MF_Y: 0d, $
	     MF_Theta: 0d, $
	     MF_Phi: 0d, $
	     axis_a: 0d, $
	     axis_b: 0d, $
	     Size: 0d, $
	     ErrSize: 0d, $
	     angle: 0d, $
	     Area: 0d, $
	     ELLIPTICITY: 0.}
if NbrObj GT 0 then begin
	     
TabObj = replicate(MyStruc, NbrTot )
Ind=0L
print, " TOTAL Number of Objects = ", NbrObj
for i=FromPatch,Nframe-1 do begin
  my_command = 'TabFace = Face'+strcompress(string(i+1), /remove_all)  
  ACK = EXECUTE( my_command) 
  for j=0,TabNbrObjPerFace[i]-1 do begin
     Obj = MyStruc
     Obj.nside = nside
     Obj.number = Ind+1
     ObjFace = TabFace[j]
     ; hs, ObjFace
  
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
     Obj.theta = Theta
     Obj.Phi = phi
     Obj.sgn = ObjFace.sgn 
     Obj.MaxVal = ObjFace.MaxVal
     Obj.Flux = ObjFace.Flux
     Obj.Size = ObjFace.Size
     if keyword_set(SZFIT) then begin
       Obj.MF_Flux = ObjFace.MF_Flux
       Obj.MF_Size = ObjFace.MF_Size
       Obj.MF_X = ObjFace.MF_X
       Obj.MF_Y = ObjFace.MF_Y
       Obj.MF_ErrFlux = ObjFace.MF_ErrFlux
       Obj.MF_ErrSize = ObjFace.MF_ErrSize
       LB = mrs_xyz2lb(PTrans, ObjFace.MF_x, ObjFace.MF_y, i)
       lb2ang, lb[0,0], lb[0,1], theta, phi
       Obj.MF_Theta = Theta
       Obj.MF_Phi = Phi
     end
     Obj.ErrorFlux = ObjFace.ErrorFlux
     Obj.Fwhm = ObjFace.Fwhm
     Obj.SNR = ObjFace.SNR
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
     Ind = Ind+1L
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
 
