;+
; NAME: 
;				PLT_COLBAR
;
; PURPOSE: 
;				draw an annotated color bar at the top of the plotting
; 			window. The window parameters are then set to the remainder of
; 			the window.
;
; CALLING:
;				plt_colbar,range,title=title,scalable=scalable,inverse=inverse,$
;				csize=csize
;
; INPUTS: 
;				range: range of the scale for the color bar
;
;OPTIONAL INPUT: 
;				title: title for the color bar
; 			scalable: use scalable pixels (useful when
;                 producing postcript files)
;       inverse: invert the color coding
;       csize: vertical color bar size (0-1)
;
; OUTPUTS:
;				an annotated color bar is drawn at the top of the window.
; 			The window parameters are then set to the remainder of the window.
; 			NOTE: to restore normal window region set !p.region=[0,0,0,0] 
;
; HISTORY:
;	Written:  A. Refregier - October 1997
;-
;-------------------------------------------------------------------------------


pro plt_colbar,range,title=title,scalable=scalable,inverse=inverse,$
  csize=csize

; set color bar vertical size
if not keyword_set(csize) then csize=.14

; reserve upper portion of the window for the color bar
!p.region=[0.,1.-csize,1.,1.]

; create a ramp for the color bar
colbar=fltarr(256,2)
colbar(*,0)=findgen(256)
colbar(*,1)=findgen(256)
if keyword_set(inverse) then colbar=-colbar
!p.charsize = 0.9	
; plot the color bar
if not keyword_set(title) then title=''
plot,[0],[0],/nodata,xstyle=12,ystyle=12
plt_image,colbar,scalable=scalable,frame=-1

plot,[0],[0],/nodata,/noerase,xrange=range,yrange=[0.,1.],$
  ystyle=12,/xstyle,title=title
oplot,[.99999,.99999]*range(0),[0.,1.]
oplot,[.99999,.99999]*range(1),[0.,1.]
!p.charsize = 1
; leave the rest of the window for the image
!p.region=[0.,0.,1.,1.-csize]

end
