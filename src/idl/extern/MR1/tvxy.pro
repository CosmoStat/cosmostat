;+
; NAME: 
;				tvxy
;
; PURPOSE: 
;				plot a vector field  map.
;
; CALLING:
; 				tvxy, DataXY, ima=ima, win=win, skip=skip, 
;                                    length_factor=length_factor, xsize=xsize, ysize=ysize,
;                                    polarization=polarization
;
; INPUT: 
; 				DataXY --- IDL 3D array [*,*,2]:  DataXY[*,*,0] = X coordinate
;                                                       [*,*,2]:  DataXY[*,*,1] = Y coordinate
;				
; KEYWORD:
;				ima --- IDL 2D array: if set, the image is first displayed, 
;                                                     and vector are overplotted
;                               skip --- scalar: the vectors are plotted for a subset of pixels
;                                                skip  (default is 10) fixes the number of pixels we skip 
;                                                between two plotted vectors
;                               length_factor --- scalar: The amplitude of the vector is plotted by
;                                                    length_factor before being plotted.
;                               win -- if set, open a new window using: window, win
;                               xsize,ysize -- scalar = size of the opened window (defaut  is 500x575)
;                               polarization -- scalar: if set, the input data corresponds a polarized 
;                                              data set where the angle between x and y axis is 45 degrees.
;                                
;
; HISTORY:
;				Written by Jean-Luc Starck, July 2007
;-
;-------------------------------------------------------------------------------

pro tvxy, xy, ima=ima, win=win, skip=skip, length_factor=length_factor, xsize=xsize, ysize=ysize, polarization=polarization, xtit=xtit, ytit=ytit, tit=tit, xtab=xtab, ytab=ytab, resol=resol

if not keyword_set(skip) then skip=10
if not keyword_set(length_factor) then length_factor=10
if not keyword_set(xsize) then  xsize=500
if not keyword_set(ysize) then  ysize=575 
if keyword_set(win) then window, win, xsize=xsize, ysize=ysize, retain=2
if not keyword_set(resol) then resol=1.

l = length_factor * resol
sz = size(xy)
n1= sz[1]
n2= sz[2]
gamma1=xy[*,*,0]
gamma2=xy[*,*,1]
if not keyword_set(resol) then resol=1
if keyword_set(ima) then k=ima $
else begin
  k = gamma1
  k[*]=0
end
help, k
plt_image, k, /frame,/col, xtit=xtit, ytit=ytit, tit=tit, xtab=xtab, ytab=ytab 
; print, resol
for j=0, n2-1, skip do begin
  for i=0, n1-1, skip do begin
    ii = float(i) * resol
    jj = float(j) * resol
    g1 = xy[i,j,0] 
    g2 = xy[i,j,1]  
    ampli = sqrt(((g1)^2) + ((g2)^2))
    if keyword_set(pol) then alph=(atan(g2,g1))/2.  $
    else alph=(atan(g2,g1))
    oplot, [ii,ii + ampli*cos(alph)*l], [jj, jj + ampli*sin(alph)*l], color=255 
    oplot, [ii,ii - ampli*cos(alph)*l], [jj, jj - ampli*sin(alph)*l], color=255 
  endfor
endfor
end


;===================

