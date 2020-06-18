;+
; NAME:
;        mrs_mr_filter
;
; PURPOSE:
;       Decompose an healpix map (nested format) into a 12 healpix face or a cube of small patches, and
;       run the MR/1 mr_filter program on each face or patch.
;
; CALLING:
;
;       mrs_mr_filter, Imag, filter, opt=opt, RmsMap=RmsMap, median=median, winsize=winsize,
;                      Patch=Patch,SizePatchDegrees=SizePatchDegrees, frac=frac
;
; INPUTS:
;     Imag -- IDL array of healpix map: Input image to be decomposed into patches 
;    
; OUTPUTS:
;     filter -- array of healpix map:Output filtered image
;                     NMaps : int = number of patchs 
;                     Nx,Ny : int = size of each patch
;                     Lon : flarr[NMaps] = Longitude if the position (Nx/2.,Ny/2.) of each map.
;                     Lat : flarr[NMaps] = Latitude if the position (Nx/2.,Ny/2.) of each map.
;                     Map : fltarr[Nx,Nx,NMaps] =  patches cube
;                     PixelSize: float = Pixel size in arc minute
;                     MapSize: float = Map size in arcmin
;                     Frac: float = overlapping factor between patches
;
; INPUT KEYWORDS:
;      Opt -- string = option to send to the mr_filter program.
;      Patch -- Scalar = if set, the image is decomposed into patches instead of of the 12 healpix faces.
;      Frac -- overlapping factor  between patches. Default is 0.05
;      SizePatchDegrees -- Size (in degrees) of each patch. Default is 10 degrees.
;      RmsMap -- IDL array of healpix map: Root mean square image.
;      median -- Scalar = if set, a simple median filtering is applied instead of the wavelet filtering
;      winsize -- Scalar = Window size of the median filtering. Default is 5.
;
; EXTERNAL CALLS:
;     mr_filter and im_filter C++ programs
;
; EXAMPLE:
;       Denoising of an image using default option
 ;             mrs_mr_filter, Imag, filter
;         
; HISTORY:
;	Written:  Jean-Luc Starck, August 2007
;--------------------------------------------------------------------------------------------------------

pro mrs_mr_filter_one_cycle, Imag, filter, opt=opt, RmsMap=RmsMap, median=median, winsize=winsize, Patch=Patch, SizePatchDegrees=SizePatchDegrees, frac=frac, ImaNoisy=ImaNoisy, FlatIma=FlatIma, msvst=msvst

if not keyword_set(Patch) then begin
  Face = H2F(Imag)
  if keyword_set(ImaNoisy) then HN = H2F(ImaNoisy)
  vs = size(imag)
  npix = vs[1]
  nside = npix2nside(npix)
  if keyword_set(RmsMap) then  RMSFace = H2F(RmsMap)
  if keyword_set(FlatIma) then  FlatFace = H2F(FlatIma)
end else begin
 PatchTrans = mrs_split(Imag, SizePatchDegrees=SizePatchDegrees, /exrec)
 Face = PatchTrans.map
 if keyword_set(RmsMap) then begin FaceRMS = mrs_split(RmsMap, SizePatchDegrees=SizePatchDegrees) & FaceRMS = FaceRMS.map & end
 if keyword_set(FlatIma) then begin  FlatFace = mrs_split(FlatIma, SizePatchDegrees=SizePatchDegrees) & FaceRMS = FaceRMS.map & end
end

if not keyword_set(opt) then opt= ' ' 
if not keyword_set(winsize) then winsize= 5
optmr= opt

optmed = ' '
if keyword_set( winsize) then optmed = '-W ' + strcompress(winsize, /remove_all)

vs= size(Face)
Nframe = vs[3]
if keyword_set(Patch) then print, "Number of patches = ", Nframe

for i=0,Nframe-1 do begin
  print, " Frame ", i+1
   f = Face[*,*,i]
   if keyword_set( median) then begin
      mr_prog, 'im_filter', f, FilterFace, opt=optmed
      ; FilterFace = median(f,  winsize) 
   end else begin
     if keyword_set(RmsMap) then begin
        RMS_f = RMSFace[*,*,i]
        namerms = gettmpfilename() + '.fits'
        writefits, namerms, RMS_f
         if not keyword_set(MSVST) then optmr= opt + ' -m 5 -R ' + namerms
     end
     if keyword_set(FlatIma) then begin
        FLAT_f = FlatFace[*,*,i]
        nameflat = gettmpfilename() + '.fits'
        writefits, nameflat, FLAT_f
        if not keyword_set(MSVST) then optmr= opt + ' -M ' + nameflat
     end
     if not keyword_set(ImaNoisy) then BEGIN
       dir = '/Users/starck/Main/BO/msvst/'
       if not keyword_set(MSVST) then mr_filter, f, FilterFace, opt=optmr $
       else mr_prog, dir+'vstestim1b ', f, FilterFace, opt=optmr   
           ;vstestimcc estimsw vstestimcc
     END else begin
       fn = HN[*,*,i]
       name_fn = gettmpfilename() + '.fits'
       writefits, name_fn, fn
       name_in = gettmpfilename() + '.fits'
       writefits, name_in, f
       optmr= opt + ' -I ' + name_in + ' -c ' + name_fn
       filename = gettmpfilename() + '.fits'
       com = "mr_pfilter " + ' ' + Optmr + ' ' +  filename
       spawn, com
       FilterFace = readfits(filename)
       delete, name_fn
       delete, name_in
       delete, filename
     end
   end
   Face[*,*,i] = FilterFace
   if keyword_set(RmsMap) then delete, namerms
end

if  not keyword_set(Patch) then  filter = f2h(Face) $
else begin
   PatchTrans.map = Face
   filter = mrs_invsplit(PatchTrans)
end
end

;===========================================


pro mrs_mr_filter, Imag, filter, opt=opt, RmsMap=RmsMap, median=median, winsize=winsize, Patch=Patch, SizePatchDegrees=SizePatchDegrees, frac=frac, ImaNoisy=ImaNoisy, FlatIma=FlatIma, i1=i1, i2=i2, i3=i3, msvst=msvst

mrs_mr_filter_one_cycle, Imag, filter, opt=opt, RmsMap=RmsMap, median=median, winsize=winsize, Patch=Patch, SizePatchDegrees=SizePatchDegrees, frac=frac, ImaNoisy=ImaNoisy, FlatIma=FlatIma, msvst=msvst

npixel = (size(Imag))[1]
nside = npix2nside(npixel)
mak_map,nside,interpole,t_interpol = 1

NCycle=4
Shift1 = [!dpi/4, -!dpi/4,  0,  0]    
Shift2 = [0, 0,           !dpi/2, -!dpi/2]
Shift3 = [0, 0,     0, 0]

TotalInterp = interpole
FilterFinal = Filter * interpole
for c=0, NCycle-1 do begin
  rotate_map_nest,interpole, Shift1[c], Shift2[c], Shift3[c], interpole2
  rotate_map_nest,Imag, Shift1[c], Shift2[c], Shift3[c], Imag2
  if keyword_set(ImaNoisy) then $
      rotate_map_nest,ImaNoisy, Shift1[c], Shift2[c], Shift3[c], ImaNoisy2
  if keyword_set(FlatIma) then $
   rotate_map_nest, FlatIma, Shift1[c], Shift2[c], Shift3[c], FlatIma2
 
 
   mrs_mr_filter_one_cycle, Imag2, filter2, opt=opt, RmsMap=RmsMap, median=median, winsize=winsize, Patch=Patch, SizePatchDegrees=SizePatchDegrees, frac=frac, ImaNoisy=ImaNoisy2, FlatIma=FlatIma2, msvst=msvst
 
   rotate_map_nest,Filter2, -Shift1[c], -Shift2[c], -Shift3[c], F2
   FilterFinal = FilterFinal + F2 * interpole2
   TotalInterp = TotalInterp + interpole2
end
   Filter  = FilterFinal / TotalInterp
end

;===========================================


;===========================================


pro test_glast, F1, F2, F3, F4
e = mrs_read('exposure_256.fits')
imag = mrs_read('counts_256.fits')
in = imag / e * 1.e10

spawn, "mr_abaque -d -e1.e-4"
opt = ' -F2 -b 3 -p -k  -n6 '
mrs_mr_filter, Imag, F1, ImaNoisy=In, opt=opt
mrs_write, "glast_fl_1e-4_F2.fits", F1
mrs_tv, F1, png='fig_glast_FL_1e-4_F2.png', title='Glast First Light', /log
f1 = rims('glast_fl_1e-4_F2.fits')

opt = ' -b 3 -p -k  -n6 '
mrs_mr_filter, Imag, F2, ImaNoisy=In, opt=opt
mrs_write, "glast_fl_1e-4.fits", F2
mrs_tv, F2, png='fig_glast_FL_1e-4.png', title='Glast First Light', /log


spawn, "mr_abaque -d -e1.e-3"
opt = ' -F2 -b 3 -p -k  -n6 '
mrs_mr_filter, Imag, F3, ImaNoisy=In, opt=opt
mrs_write, "glast_fl_1e-3_F2.fits", F3
mrs_tv, F3, png='fig_glast_FL_1e-3_F2.png', title='Glast First Light', /log

opt = ' -b 3 -p -k  -n6 '
mrs_mr_filter, Imag, F4, ImaNoisy=In, opt=opt
mrs_write, "glast_fl_1e-3.fits", F4
mrs_tv, F4, png='fig_glast_FL_1e-3.png', title='Glast First Light', /log

end



;===========================================

pro test_glast1, F1, F2, F3
e = mrs_read('expo_100_300_256.fits')
imag = mrs_read('counts_100_300_256.fits')
in = imag / e * 1.e10

spawn, "mr_abaque -d -e1.e-6"
opt = ' -F2 -b 3 -p -k -o -n6 '
mrs_mr_filter, Imag, F1, ImaNoisy=In, opt=opt
mrs_write, "glast_fl_100-300_1e-6_F2.fits", F1
mrs_tv, F1, png='fig_glast_100-300_FL_1e-6_F2.png', title='Glast First Light (100-300)', /log

e = mrs_read('expo_300_1000_256.fits')
imag = mrs_read('counts_300_1000_256.fits')
in = imag / e * 1.e10
mrs_mr_filter, Imag, F2, ImaNoisy=In, opt=opt
mrs_write, "glast_fl_300-1000_1e-6_F2.fits", F2
mrs_tv, F2, png='fig_glast_300-1000_FL_1e-6_F2.png', title='Glast First Light (300-1000)', /log

e = mrs_read('expo_sup1000_256.fits')
imag = mrs_read('counts_sup1000_256.fits')
in = imag / e * 1.e10
mrs_mr_filter, Imag, F3, ImaNoisy=In, opt=opt
mrs_write, "glast_fl_sup1000_1e-6_F2.fits", F3
mrs_tv, F3+1, png='fig_glast_sup1000_FL_1e-6_F2.png', title='Glast First Light (sup1000)', /log


end



;===========================================

pro test_glast2, F1, F2, F3


spawn, "mr_abaque -d -e1.e-6"
opt = ' -f3 -F2 -b 3 -p -k  -n6 '

e = mrs_read('expo_sup1000_256.fits')
imag = mrs_read('counts_sup1000_256.fits')
in = imag / e * 1.e10
mrs_mr_filter, Imag, F1, ImaNoisy=In, opt=opt
mrs_write, "glast_fl_sup1000_1e-6_F2.fits", F1
mrs_tv, F1, png='fig_glast_sup1000_FL_1e-6_F2.png', title='Glast First Light (sup1000)', /log


e = mrs_read('expo_sup1000_512.fits')
imag = mrs_read('counts_sup1000_512.fits')
in = imag / e * 1.e10
mrs_mr_filter, Imag, F2, ImaNoisy=In, opt=opt
mrs_write, "glast_fl_sup1000_512_1e-6_F2.fits", F2
mrs_tv, F2, png='fig_glast_sup1000_512_FL_1e-6_F2.png', title='Glast First Light (nside=512, [1GeV,100GeV])', /log

end


;===========================================

function add_wtresi, d, r, ns1=ns1, ns2=ns2, wr=wr, coef=coef, sr=sr
d1  = dblarr(12L*256L^2)
d1[*] = min(d)
ind = where(d GT 0 and d LT 10000, c)
if c GT 0 then d1[ind] = d[ind]
if not keyword_set(ns1) then ns1=0
if not keyword_set(ns2) then ns2=5
if not keyword_set(coef) then coef=0.

if not keyword_set(wr) then mrs_wttrans, r, wr, nbrscale=7
sr = wr.coef[*,ns1:ns2]
f = wr.coef[*,ns1:ns2]*coef
return, d+f
end

;===========================================

pro plot_add_resi

e = mrs_read('exposure_256.fits')
imag = mrs_read('counts_256.fits')
in = imag / e * 1.e10

; spawn, "mr_abaque -d -e1.e-4"
; opt = ' -F2 -b 3 -p -k  -n6 '
; mrs_mr_filter, Imag, F1, ImaNoisy=In, opt=opt
; mrs_write, "glast_fl_1e-4_F2.fits", F1
; mrs_tv, F1, png='fig_glast_FL_1e-4_F2.png', title='Glast First Light', /log
f1 = rims('glast_fl_1e-4_F2.fits')
d1  = dblarr(12L*256L^2)
d1[*] = min(f1)
ind = where(f1 GT 0 and f1 LT 10000, c)
if c GT 0 then d1[ind] = f1[ind]
r = in - d1
wr=0
tvs, add_wtresi(d1, r, wr=wr,coef=0.1), /log
tvs, d1 + r*0.1, png='large_band_resi_0p1.png', title='Large band, WT denoising + 0.1 Residual', /log
tvs, d1 + r*0.2, png='large_band_resi_0p2.png', title='Large band, WT denoising + 0.2 Residual', /log
mrs_write, 'glast_fl_1e-4_F2_resi_0p1.fits', d1 + r*0.1
mrs_write, 'glast_fl_1e-4_F2_resi_0p2.fits', d1 + r*0.2
; tvs, d1 + r*0.5, png='large_band_resi_0p5.png', title='Large band, WT denoising + 0.5 Residual', /log


; RESI 100 - 300
e = mrs_read('expo_100_300_256.fits')
imag = mrs_read('counts_100_300_256.fits')
in = imag / e * 1.e10
; spawn, "mr_abaque -d -e1.e-6"
; opt = ' -F2 -b 3 -p -k -o -n6 '
; mrs_mr_filter, Imag, F1, ImaNoisy=In, opt=opt
F1 = mrs_read("glast_fl_100-300_1e-6_F2.fits")
d1  = dblarr(12L*256L^2)
d1[*] = min(f1)
ind = where(f1 GT 0 and f1 LT 10000, c)
if c GT 0 then d1[ind] = f1[ind]
r = in - d1
wr=0
tvs, add_wtresi(d1, r, wr=wr,coef=0.1), /log
mrs_write, "glast_fl_100-300_1e-6_F2_resi_0p1.fits",  d1 + r*0.1
mrs_tv,  d1 + r*0.1, png='fig_glast_100-300_FL_1e-6_F2_resi_0p1.png', title='Glast First Light (100-300) + 0.1 Residual', /log
mrs_write, "glast_fl_100-300_1e-6_F2_resi_0p2.fits",  d1 + r*0.2
mrs_tv,  d1 + r*0.2, png='fig_glast_100-300_FL_1e-6_F2_resi_0p2.png', title='Glast First Light (100-300) + 0.2 Residual', /log


; RESI 300 - 1000
e = mrs_read('expo_300_1000_256.fits')
imag = mrs_read('counts_300_1000_256.fits')
in = imag / e * 1.e10
; mrs_mr_filter, Imag, F2, ImaNoisy=In, opt=opt
; mrs_write, "glast_fl_300-1000_1e-6_F2.fits", F2
F1 = mrs_read("glast_fl_300-1000_1e-6_F2.fits")
mrs_tv, F1, png='fig_glast_300-1000_FL_1e-6_F2.png', title='Glast First Light (300-1000)', /log
d1  = dblarr(12L*256L^2)
d1[*] = min(f1)
ind = where(f1 GT 0 and f1 LT 10000, c)
if c GT 0 then d1[ind] = f1[ind]
r = in - d1
wr=0
tvs, add_wtresi(d1, r, wr=wr,coef=0.01), /log
mrs_tv,  d1 + r*0.1, png='fig_glast_300-1000_FL_1e-6_F2_resi_0p1.png', title='Glast First Light (300-1000) + 0.1 Residual', /log
mrs_write, "glast_fl_300-1000_1e-6_F2_resi_0p1.fits", d1 + r*0.1
mrs_tv,  d1 + r*0.2, png='fig_glast_300-1000_FL_1e-6_F2_resi_0p2.png', title='Glast First Light (300-1000) + 0.2 Residual', /log
mrs_write, "glast_fl_300-1000_1e-6_F2_resi_0p2.fits", d1 + r*0.2

; RESI SUP 1000
e = mrs_read('expo_sup1000_256.fits')
imag = mrs_read('counts_sup1000_256.fits')
in = imag / e * 1.e10
; mrs_mr_filter, Imag, F3, ImaNoisy=In, opt=opt
; mrs_write, "glast_fl_sup1000_1e-6_F2.fits", F3
; mrs_tv, F3+1, png='fig_glast_sup1000_FL_1e-6_F2.png', title='Glast First Light (sup1000)', /log
F1 = mrs_read("glast_fl_sup1000_1e-6_F2.fits")
d1  = dblarr(12L*256L^2)
d1[*] = min(f1)
ind = where(f1 GT 0 and f1 LT 10000, c)
if c GT 0 then d1[ind] = f1[ind]
r = in - d1
wr=0
tvs, add_wtresi(d1, r, wr=wr,coef=0.1), /log
r1 = add_wtresi(d1, r, wr=wr,coef=0.1, ns1=1, sr=sr)
mrs_tv,  d1 + sr*0.1, png='glast_fl_sup1000_1e-6_F2_resi_0p1G.png', title='Glast First Light (sup1000) + 0.1 Residual * G ', /log
mrs_write, "glast_fl_sup1000_1e-6_F2_resi_0p1G.fits", d1 + r*0.1
mrs_tv,  d1 + sr*0.2, png='glast_fl_sup1000_1e-6_F2_resi_0p2G.png', title='Glast First Light (sup1000) + 0.2 Residual * G ', /log
mrs_write, "glast_fl_sup1000_1e-6_F2_resi_0p2G.fits", d1 + r*0.2

mrs_tv,  d1 + r*0.1, png='glast_fl_sup1000_1e-6_F2_resi_0p1.png', title='Glast First Light (sup1000) + 0.1 Residual', /log
mrs_write, "glast_fl_sup1000_1e-6_F2_resi_0p1.fits", d1 + r*0.1
mrs_tv,  d1 + r*0.2, png='glast_fl_sup1000_1e-6_F2_resi_0p2.png', title='Glast First Light (sup1000) + 0.2 Residual', /log
mrs_write, "glast_fl_sup1000_1e-6_F2_resi_0p2.fits", d1 + r*0.2


end

;===========================================

pro glast_filter, i, e, Filter, wresi=wresi, resi=resi, InitF=InitF, opt=opt, BgrResi= BgrResi
if not keyword_set(opt) then opt=' -v -f3 -A -m10 -E1e-6 -v -p -k  -n7 -K'

mrs_mr_filter, i, F1,  opt=opt, i1=i1, i2=i2, i3=i3
InitF = F1

resi = i - F1
mrs_wttrans, resi, w, nbrscale=7
wresi=w
bgr = w.coef[*,6]+w.coef[*,5]+w.coef[*,4]

bgr2 = total( w.coef[*,1:3], 2)

F2 = F1 + bgr

Filter=F2
end


;===========================================

pro test_glast3, F1, i1=i1, i2=i2, i3=i3


opt = ' -v -f3 -A -m10 -E1e-6 -v -p -k  -n7 -K'
; opt='  -v -s4.7 -M0 -n7 -p'

e = mrs_read('expo_sup1000_256.fits')
imag = mrs_read('counts_sup1000_256.fits')
Flat = 1.e10 / e 
; mrs_mr_filter, Imag, F1, FlatIma=Flat, opt=opt,i1=i1, i2=i2, i3=i3 

mrs_mr_filter, Imag, F1,  opt=opt, i1=i1, i2=i2, i3=i3

resi = imag - F1
mrs_wttrans, resi, w, nbrscale=7
bgr = w.coef[*,6]+w.coef[*,5]+w.coef[*,4]

bgr2 = total( w.coef[*,1:3], 2)

F2 = F1 + bgr
mrs_write, "glast_new_fl_sup1000_1e-6_F2.fits", F2 /flat
mrs_tv, F2/flat, png='fig_glast_new_sup1000_FL_1e-6_F2.png', title='Glast First Light (sup1000)', /log


end


;===========================================

pro test_face, i=i, RF

e = mrs_read('expo_sup1000_256.fits')
imag = mrs_read('counts_sup1000_256.fits')
Flat = 1.e10 / e 
hf = h2f(Flat)
hi = h2f(imag)
rf = hi
for i=0, 11 do begin
flat1 = hf[*,*,i]
i1 = hi[*,*,i]
writefits, 'tmp_flat.fits', flat1
opt = ' -f3 -F2 -A -m10 -E1e-6 -v -p -k  -n6 -M tmp_flat.fits -G0.1 '
mr_filter, i1, r1, opt=opt
wset, 0
load, i1
wset, 1
load, r1
info, r1
RF[*,*,i] = r1
end

end


