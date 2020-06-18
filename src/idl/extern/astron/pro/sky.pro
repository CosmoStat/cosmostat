pro sky,image,skymode,skysig,CIRCLERAD=circlerad,SILENT=silent
;+
; NAME:
;       SKY
; PURPOSE:
;       Determine the sky level in an image using the the procedure MMM
; EXPLANATION:
;       Approximately 4000 uniformly spaced pixels are selected for the
;       computation.  Adapted from the DAOPHOT routine of the same name.
;
; CALLING SEQUENCE:
;       SKY, image, [ skymode, skysig, CIRCLE = ,/SILENT ]
;
; INPUTS:
;       IMAGE - One or two dimensional array
;
; OPTIONAL OUTPUT ARRAYS:
;       SKYMODE - Scalar, giving the mode of the sky pixel values of the 
;               array IMAGE, as determined by the procedure MMM.
;       SKYSIG -  Scalar, giving standard deviation of sky brightness
;
; INPUT KEYWORD PARAMETERS:
;       CIRCLERAD - Use the keyword to have SKY only select pixels within
;               the specified pixel radius of the center of the image.  If 
;               CIRCLERAD = 1, then the radius is set equal to half the image
;               width.   Useful when the data is restricted to a circular area
;               of the image.
;
;       /SILENT - If this keyword is supplied and non-zero, then SKY will not
;               display the sky value and sigma at the terminal
;
; PROCEDURE:
;       A regular grid of points, not exceeding 4000 in number, is extracted
;       from the image array.  The mode of these pixel values is determined
;       by the procedure MMM.
;
; PROCEDURE CALLS:
;       MMM
; REVISION HISTORY:
;       Written, W. Landsman   STX Co.            September, 1987     
;       Changed INDGEN to LINDGEN                 January, 1994
;       Fixed display of # of points used         March, 1994
;       Converted to IDL V5.0   W. Landsman   September 1997
;-
 On_error,2              ;Return to caller
 maxsky = 4000          ;Maximum # of pixels to be used in sky calculation

 if N_params() eq 0 then begin
        print,'Syntax - sky, image, [ skymode, skysig , CIRCLE = , /SILENT ]'
        return
 endif

 s = size(image)                                  
 case s[0] of                           ;Number of dimensions in array? 

 1: begin                                       ; 1-d vector
    istep = long(s[1]/float(maxsky) + 0.5) > 1
    npts = s[1]/istep 
    ilow = (s[1] - istep*npts)/2
    index = ilow + lindgen(npts)*istep           ;Modified Jan. 94
    end

 2: begin                                       ;2-D array
    if keyword_set(CIRCLERAD) then begin 
       width = max([s[1],s[2]])
       if circlerad EQ 1 then rad = width/2 else rad = long(circlerad)
       npts = !PI*rad^2
    endif else npts = s[4]
    fractn = float(maxsky)/npts
    istep = fix(0.5 + 1/sqrt(fractn)) > 1
    ilow = (s[1] - istep*((s[1] - 2)/istep) -1) / 2
    n_x = s[1] / istep
    jstep = s[2] *n_x/maxsky + 1
    jlow = (s[2] - jstep*((s[2] -2)/jstep) - 1) /2
    n_y = s[2]/jstep
    npts = n_x*n_y
    yy  =  jlow + jstep*(lindgen(npts)/n_x)
    xx =  ilow + istep*(lindgen(npts) mod n_x)
    index = xx + s[1]*yy
    if keyword_set(CIRCLERAD) then begin
      pixrad = (xx-s[1]/2)^2 + (yy-s[2]/2)^2
      good = where(pixrad lt rad^2,ngood)
      if ngood GT 0 then index = index[good] else $
          message,'ERROR - Too few pixels within specified circle radius'
    endif
   end 

 else: message,'Input array (first parameter) must be 1 or 2 dimensional'

 endcase

 mmm, image[index], skymode, skysig

 skymode = float(skymode)  &  skysig = float(skysig)
 if not keyword_set(SILENT) then begin
        print,'Number of points used to find sky = ',N_elements(index)
        print,'Approximate sky value for this frame = ',skymode
        print,'Standard deviation of sky brightness = ',skysig
 endif

 return
 end
