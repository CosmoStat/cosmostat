;+
; NAME: 
;				MK_JPEG
;
; PURPOSE:
;				Produce a jpeg image by capturing the current window.
; 				Alternatively, an image array can be turned into a jpeg image.
; 			The current color table is used
;	
;
; CALLING:
;				mk_jpeg,image=m,fname=fname,quality=quality,progressive=progressive,$
;       inframe=inframe
;
; INPUTS:
;				current window unless the image parameter is used
;       current color table
; 
; OPTIONAL INPUT: 
;				image: image array; if not provided, the current
;                        window is captured. 
;				fname: output file name. (default='idl.jpeg')
;       quality: jpeg quality factor (0: bad, 100:excellent,
;                default=75)
;       progressive: jpeg progressice storage (image will
;                    reconstruct progressively when loaded)
;      inframe: capture only the portion of the window
;               which is inside the axis frame. This is
;               useful to capture plots-image overlays
;               produced by plt_image, while avoiding spaces
;               on the sides.	
;
; OUTPUTS: 
;				window (or image) saved as a true-color jpeg image
; 			capture the current window if an image is not provided
;
; HISTORY:
;	August 5, 1997 - Written by A. Refregier
;-
;-------------------------------------------------------------------------------


pro mk_jpeg,image=m,fname=fname,quality=quality,progressive=progressive,$
            inframe=inframe, inwin=inwin

if keyword_set(inwin) then begin
   write_jpeg, fname, tvrd(true=1), true=1
end else BEGIN

; capture the current window if an image is not provided
if not keyword_set(m) then begin
  m=tvrd(0) 
  if keyword_set(inframe) then begin    ; capture portion inside frame only
     px=!x.window*!d.x_vsize     ; frame limits
     py=!y.window*!d.y_vsize 
     m=m(px(0):px(1),py(0):py(1))
  endif
endif

; compute image dimensions and characteristics
n_m1=n_elements(m(*,0))
n_m2=n_elements(m(0,*))
m_max=float(max(m)) & m_min=float(min(m))
m_ran=m_max-m_min

; load in current color tables
tvlct,rc,gc,bc,/get
n_rc=n_elements(rc)
n_gc=n_elements(gc)
n_bc=n_elements(bc)

; compute the three color arrays by interpolating the color tables
r=rc  (fix(((float(m)-m_min)/m_ran)*n_rc))
b=bc  (fix(((float(m)-m_min)/m_ran)*n_bc))
g=gc  (fix(((float(m)-m_min)/m_ran)*n_gc))

; write to jpeg file as a true color image
if not keyword_set(fname) then fname='idl.jpeg'
write_jpeg,fname,[[[r]],[[g]],[[b]]],true=3,$
  quality=quality,progressive=progressive
END

end

