;+
; NAME: 
;				mrsp_tvf
;
; PURPOSE: 
;				plot one  Healpix face of a polarized vector field  map 
;
; CALLING:
; 				mrsp_tvf, map, face, win=win, skip=skip, length_factor=length_factor, xsize=xsize, ysize=ysize,
;                                        log=log, maxval=maxval, minval=minval, xtit=xtit, ytit=ytit, tit=tit, jpg=jpg, deg=deg
;
; INPUT: 
; 				map --- IDL 2D array Healpix POLARIZED map in NESTED format:  map[*,0] = T field,	map[*,1] = Q field,		map[*,2] = U field
;				face --- int: face number between 0 and 11 
;				
; KEYWORD:
;				ima --- IDL 2D array: if set, the image is first displayed, 
;                                                     and vector are overplotted
;                               skip --- scalar: the vectors are plotted for a subset of pixels
;                                                skip  (default is 10) fixes the number of pixels we skip 
;                                                between two plotted vectors
;                               length_factor --- scalar: The amplitude of the vector is plotted by
;                                                    length_factor before being plotted.
;                               win -- if set, open a new window using: window, win, xsize=xsize, ysize=ysize
;                               xsize,ysize -- scalar = size of the opened window (defaut  is 500x575)
;                               xtit -- string: x-axis title
;                               ytit -- string: y-axis title
;                               tit -- string:  title of the plot
;                               jpg -- string: if set, a jpeg file is created with filename given by jpg
;                               log -- scalar: if set, plot the temperature in log
;                               minval, maxval -- scalar: min value and max value of the temperature map to be plotted
;                               deg -- scalar: if set, x and y axis are in degree (arcminute by default)
;                                
;
; HISTORY:
;				Written by Jean-Luc Starck, July 2007
;-
;-------------------------------------------------------------------------------

pro mrsp_tvf, map, face, win=win, skip=skip, length_factor=length_factor, xsize=xsize, ysize=ysize, log=log, maxval=maxval, minval=minval, xtit=xtit, ytit=ytit, tit=tit, jpg=jpg, deg=deg

if not keyword_set(skip) then skip=10
if not keyword_set(length_factor) then length_factor=10
if not keyword_set(xsize) then  xsize=500
if not keyword_set(ysize) then  ysize=575 

nside= gnside(map[*,0])

S1 = H2F(map[*,0])
S2 = H2F(map[*,1])
S3 = H2F(map[*,2])
sz = size(S1)
n1= sz[1]
n2= sz[2]
temp = S1[*,*,face]
if keyword_set(maxval) then temp = temp < maxval
if keyword_set(minval) then temp = temp > minval

if keyword_set(log) then temp = alog(temp+1)

Mxy = fltarr(n1,n2, 2)
Mxy[*,*,0] = S2[*,*,face]
Mxy[*,*,1] = S3[*,*,face]
; help, Mxy
FigTit = 'Healpix face ' + strc(face+1)
if keyword_set(tit) then FigTit= FigTit + ':  ' + tit
if keyword_set(log) then FigTit= 'LOG ' + FigTit 

resol = pixel_size(nside)
xtit='arc minute'
ytit='arc minute'
xtab = findgen(n1)  * resol
ytab = findgen(n2)  * resol
if keyword_set(deg) then begin
  xtit='degrees'
  ytit='degrees'
  xtab = xtab / 60.
  ytab = ytab / 60.
  resol = resol / 60.
end


tvxy, Mxy, ima=temp, win=win, skip=skip, length_factor=length_factor, xsize=xsize, ysize=ysize, /polarization, xtit=xtit, ytit=ytit, tit=FigTit, xtab=xtab, ytab=ytab, resol=resol

if keyword_set(jpg) then mk_jpeg, quality=100, fname=jpg

end
;   mrsp_tvf, ip, 5, skip=5, len=30, /log, tit='test', win=1, /deg

;===================

pro mrp_tvf, Ima, win=win, skip=skip, length_factor=length_factor, xsize=xsize, ysize=ysize, log=log, maxval=maxval, minval=minval, xtit=xtit, ytit=ytit, tit=tit, jpg=jpg, deg=deg

sz = size(Ima)
n1= sz[1]
n2= sz[2]

if not keyword_set(skip) then skip=10
if not keyword_set(length_factor) then length_factor=10
if not keyword_set(xsize) then  xsize=500
if not keyword_set(ysize) then  ysize=575 

nside=n1

temp = Ima[*,*,0]
if keyword_set(maxval) then temp = temp < maxval
if keyword_set(minval) then temp = temp > minval

if keyword_set(log) then temp = alog(temp+1)
 
Mxy = Ima[*,*,1:2]
; help, Mxy
FigTit = ' '  
if keyword_set(tit) then FigTit= FigTit + ':  ' + tit
if keyword_set(log) then FigTit= 'LOG ' + FigTit 

resol = pixel_size(nside)
xtit='arc minute'
ytit='arc minute'
xtab = findgen(n1)  * resol
ytab = findgen(n2)  * resol
if keyword_set(deg) then begin
  xtit='degrees'
  ytit='degrees'
  xtab = xtab / 60.
  ytab = ytab / 60.
  resol = resol / 60.
end

tvxy, Mxy, ima=temp, win=win, skip=skip, length_factor=length_factor, xsize=xsize, ysize=ysize, /polarization, xtit=xtit, ytit=ytit, tit=FigTit, xtab=xtab, ytab=ytab, resol=resol

if keyword_set(jpg) then mk_jpeg, quality=100, fname=jpg

end

;===================

pro ppp, ip, win=win
mrsp_tvf, ip, 5, skip=5, len=30, /log, tit='test', win=win, /deg, win=1
end
;% RESTORE: Restored variable: DUST.
;% RESTORE: Restored variable: DUST_3.
;% RESTORE: Restored variable: DUST_5.
;% RESTORE: Restored variable: DUST_TEB3.
;% RESTORE: Restored variable: DUST_TEB5.
;% RESTORE: Restored variable: DUST_TQU3.
;% RESTORE: Restored variable: DUST_TQU5.
