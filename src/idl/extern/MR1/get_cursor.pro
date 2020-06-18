pro get_cursor, hd, pict, x0=x0, y0=y0, Filename=Filename, zoom=zoom
;+
; NAME:
;	GET_CURSOR
; PURPOSE:   
;	Display values associated with the cursor.
;	GET_CURSOR displays different information depending whether the user 
;	supplied an image array, and/or a FITS header array.
;       This routine is extracted from the Astronomical User Library, and has
;       been modified. Values are printed only when the user click on the 
;       left mouse button, and a zoom factor has been introduced for the
;       case where the image in the window has been zoomed.
;	X and Y values, etc., of the pixel under the cursor are constantly
;	displayed.  
;	Pressing left or center mouse button prints a line of output, and 
;	starts a new line.
;	Pressing right mouse button exits the procedure.
;	If the keyword FILENAME is defined, the date and time, and a heading 
;	will be printed in the file before the data.
;;
;
; CALLING SEQUENCE(S):
;	GET_CURSOR          ;Display x,y and byte intensity (inten)
;	
;	GET_CURSOR, pict   ;Display x,y,inten, and also pixel value (from image array)
;	
;	GET_CURSOR, hdr, pict, [X0 =, Y0 =, FILE =]        
;
;
; OPTIONAL INPUTS:
;	Hdr  -- IDL structure : FITS Header array
;	pict  -- array : Array containing values that are displayed.  Any type.
;
; OPTIONAL KEYWORD INPUTS:
;	X0,Y0  -- scalars :  coordinates of lower left-hand corner of input image 
;		If not supplied then GET_CURSOR assumes that pixel (0,0) on the
;		image is located at position (0,0) in the window 
;	FILENAME  -- string:  name of file to where GET_CURSOR data can be saved.
;		Data will only be saved if left or center mouse button
;		are pressed.
;
;        ZOOM -- scalar :  Zoom facted between the display image and input image
;
; EXAMPLE:
;       GET_CURSOR          ;Display x,y and byte intensity (inten)
;
;       GET_CURSOR, pict   ;Display x,y,inten, and also pixel value (from image array)
;
;       GET_CURSOR, hdr, pict, [X0 =, Y0 =, FILE =]
;
; MODIFICATION HISTORY:
;	Written,  K. Rhode,  STX  May 1990
;	Added keyword FILENAME	D. Alexander  June 1991
;       Don't write to Journal file   W. Landsman    March 1993
;	Use astrometry structure  W. Landsman      Feb 1994
;  Modified for Mac IDL          I.   Freedman     April 1994
;   13-MAR-1996 J. SAM LONE: update header
;   13-OCT-1998 JL Starck : bracket instead of parenthesis
;
;-
 
 On_error,2    ;if an error occurs, return to caller
if not keyword_set(zoom) then zoom = 1
all = 0

 if !VERSION.OS EQ 'MacOS' then begin
  print,'Macintosh Mouse maps to "LEFT" button and LEFT and RIGHT arrow keys to'
  print,'"CENTER" and "RIGHT" buttons respectively'
 endif

 f_header = 0b           ;True if a FITS header supplied
 PictPar =  0b           ;True if an image array supplied
 f_astrom = 0b           ;True if FITS header contains astrometry
 f_bscale = 0b           ;True if FITS header contains BSCALE factors
 f_imhd   = 0b           ;True if image array is in HD (1 parameter)
 npar = N_params()
 fileflag=0		;True once left or middle mouse button pressed

 if !D.WINDOW EQ -1 then begin
	message,'ERROR - No image window active',/INF
	return
 endif

 if !VERSION.OS  NE 'MacOS' then begin
  openw,term,filepath(/TERMINAL),/GET_LUN,/STREAM         ;Output terminal
 endif else  term = -1
;
; Determine if there is an X or Y offset in common block or supplied keyword.
;
if not keyword_set(X0) then x0 = 0
if not keyword_set(Y0) then y0 = 0
;
if (!D.FLAGS and 256) EQ 256 then wshow,!D.WINDOW  ;Bring active window to foreground
;
; Print formats and header for different astrometry,image, BSCALE combinations
;
cr = string("15b)
line0 = '  X     Y     Byte Inten'
line1 = '  X     Y     Byte Inten   Value'
line2 = '  X     Y     Byte Inten        RA                Dec   '
line3 = '  X     Y   ByteInten   Value        RA             Dec         Flux' 
line4 = '  X     Y   ByteInten   Value        RA             Dec' 
line5 = '  X     Y   ByteInten   Value   Flux'
f0 = "($,i4,2x,i4,6x,i4,a)"
f1 = "($,i4,2x,i4,6x,i4,3x,a,a)"
f2 = "($,i4,2x,i4,6x,i4,5x,i4,i4,1x,f6.2,3x,i4,i4,1x,f6.2,a)"
f3 = "($,i4,2x,i4,2x,i4,3x,a,2x,i4,i4,1x,f6.2,2x,i4,i4,1x,f6.2,5x,e8.2,a)"
f4 = "($,i4,2x,i4,2x,i4,3x,a,2x,i4,i4,1x,f6.2,2x,i4,i4,1x,f6.2,a)"
f5 = "($,i4,2x,i4,2x,i4,3x,a,5x,e8.2,a)"

if (npar gt 0) then begin
  type = size(hd)
  if (npar eq 1) and (type[0] eq 2) then begin
    PictPar = 1b  & f_imhd = 1b 
    imtype = type
    Ny = imtype[2]
  endif else if (type[2] ne 7) or (type[0] ne 1) then begin
    print,'Calling sequence options: GETCURSOR,'
    print,'                          GETCURSOR,IM  where IM is a 2-D image,'
    print,'                          GETCURSOR,HD  where HD is a FITS header,'
    print,'                          or  GETCURSOR,HD,IM.'
    return
  endif else if (type[2] eq 7) and (type[0] eq 1) then f_header = 1b
  if (npar eq 2) then begin
    PictPar = 1b & f_header = 1b
    imtype = size(pict)
    Ny = imtype[2]
    if (imtype[0] lt 2) or $
     (imtype[imtype[0]+2] ne imtype[1]*imtype[2]) then $
       message,'Image array (second parameter) is not two dimensional.'
  endif
endif    


;
; get information from the header
;
if f_header then begin     

  extast,hd,astr,noparams
  if (noparams ge 0) then f_astrom = 1b

  if PictPar then begin
  bscale = sxpar(hd,'BSCALE')
  if (bscale ne 0) then begin
    bzero = sxpar(hd,'BZERO')
    bunit = sxpar(hd,'BUNIT')
    if !ERR ge 0 then $ 
    if f_astrom then line3 = line3 + '('+bunit+ ')' else $
                     line5 = line5 + '('+bunit+')'
    f_bscale = 1b
  endif
  endif
endif
;
print,'Press left or center mouse button for new output line,'
print,'... right mouse button to exit.  
;
; different print statements, depending on the parameters
case 1 of

(PictPar eq 0b) and (f_astrom eq 0b):  begin   
   curtype = 0 & print, line0  & end      ;No image or header info

(PictPar) and (f_astrom eq 0b) and (f_bscale eq 0b): begin
   curtype = 1  & print,line1 & end       ;Only image array supplied

(PictPar eq 0b) and (f_astrom) and (f_bscale eq 0b): begin 
   curtype = 2  & print,line2 & end       ;Astrometry but no image array

(PictPar) and (f_astrom) and (f_bscale): begin
   curtype =3   & print,line3 & end       ;Image array + astrometry + BSCALE

(PictPar) and (f_astrom) and (f_bscale eq 0b): begin
   curtype = 4  & print,line4 & end       ;Image array +astrometry

(PictPar) and (f_astrom eq 0b) and (f_bscale): begin
   curtype = 5  & print,line5 & end       ;Image array + BSCALE

endcase
;
;
LOOP: sv_err = !ERR
!ERR = 0

cursor,x,y,2,/DEVICE,/CHANGE                                 
 if !ERR EQ 4 then begin
    if !VERSION.OS NE 'MacOS' then free_lun,term
   
    if fileflag AND !VERSION.OS NE 'MacOS' then free_lun,lun
    return
 endif

 if (!ERR GT 1) and (!ERR NE sv_err) then begin  ; have left or center buttons been pressed?
    print,form="($,a)",string("12b)   ; print a form feed
    if keyword_set(filename) and (not fileflag) then begin	
                                 ; open file & print table header to file
	get_lun,lun
	openw,lun,filename
	printf,lun,'GETCURSOR:   ',systime(0)	;print time and date to file
        case 1 of  		;different print statements for file, depending on parameters

        (PictPar eq 0b) and (f_astrom eq 0b) : begin
           printf, lun, line0  & end			;No image or header info

	(PictPar) and (f_astrom eq 0b) and (f_bscale eq 0b) : begin
	   printf, lun, line1 & end			;Only image array supplied

	(PictPar eq 0b) and (f_astrom) and (f_bscale eq 0b) : begin
   	   printf, lun, line2 & end			;Astrometry but no image array

	(PictPar) and (f_astrom) and (f_bscale) : begin
    	   printf, lun, line3 & end			;Image array + astrometry + BSCALE

	(PictPar) and (f_astrom) and (f_bscale eq 0b) : begin
   	   printf, lun, line4 & end			;Image array + astrometry

	(PictPar) and (f_astrom eq 0b) and (f_bscale) : begin
   	   printf, lun, line5 & end			;Image array + BSCALE
	endcase
 	fileflag=1
    endif
    if keyword_set(filename) then begin
	case curtype of 
	   0: printf, lun, form=f0, x, y, inten
	   1: printf, lun, form=f1, x, y, inten, value 
	   2: printf, lun, form=f2, x, y, inten, ihr, imin, xsec, deg, mn, sc
	   3: printf, lun, form=f3, x, y, inten, value, ihr, imin, xsec, deg, mn, sc, flux
	   4: printf, lun, form=f4, x, y, inten, value, ihr, imin, xsec, deg, mn, sc
	   5: printf, lun, form=f5, x, y, inten, value, flux
	endcase
    endif
 endif

if (!ERR EQ 1) OR (all EQ 1) then BEGIN
  x = x>0 & y = y>0
  inten = fix(tvrd(x,y,1,1))   ; read the byte intensity 
if !order EQ 1 then y = Ny*zoom - 1 - y
if (x - x0) GE 0 then x = fix((x - x0) / zoom) $
else x = fix((x - x0) / zoom-1)

if (y - y0) GE 0 then y = fix((y - y0) / zoom) $
else y = fix((y - y0) / zoom-1)

 if PictPar EQ 1b then begin
      if (x LT 0) or (x GE imtype[1]) or $
         (y LT 0) or (y GE imtype[2]) then value = 0 else $
      if f_imhd then value = hd[x,y] else value = pict[x,y]
 endif

 if f_astrom then begin
      xy2ad, x, y, astr, a, d                      ; convert to ra and dec
      radec, a, d, ihr, imin, xsec, deg, mn, sc   ; change to hr,min format
 endif

 if f_bscale  then flux = bscale*value + bzero  

 case curtype of
	0:  printf,term,form=f0,x,y,inten,cr  
	1:  printf,term,form=f1,x,y,inten,value,cr 
	2:  printf,term,form=f2,x,y,inten,ihr,imin,xsec,deg,mn,sc,cr        
	3:  printf,term,form=f3,x,y,inten,value,ihr,imin,xsec,deg,mn,sc,flux,cr
	4:  printf,term,form=f4,x,y,inten,value,ihr,imin,xsec,deg,mn,sc,cr 
	5:  printf,term,form=f5,x,y,inten,value,flux,cr
 endcase
END
 goto,LOOP

 end
