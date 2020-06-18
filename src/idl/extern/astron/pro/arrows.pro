pro arrows,h,xcen,ycen,thick=thick,charsize=charsize,arrowlen=arrowlen, $
             color=color,NotVertex=NotVertex,Normal = normal,Data=data,font=font
;+
; NAME:
;      ARROWS
; PURPOSE:
;      To display "weathervane" directional arrows on an astronomical image 
; EXPLANATION:
;      Overlays a graphic showing orientation of North and East.
;
; CALLING SEQUENCE:
;      ARROWS,h, [ xcen, ycen, ARROWLEN= , CHARSIZE=  COLOR= , /DATA
;                              FONT=, /NORMAL, /NOTVERTEX, THICK=  ]
;
; INPUTS:
;       h - FITS or STSDAS header array, must include astrometry
;
; OPTIONAL INPUTS:
;       xcen,ycen - numeric scalars, specifying the center position of
;		arrows.   Position in device units unless the /NORMALIZED 
;		keyword is specified.   If not supplied, then ARROWS
;		will prompt for xcen and ycen
;
; OPTIONAL KEYWORD INPUTS:
;       arrowlen  - length of arrows in terms of normal Y size of vector-drawn
;                     character,  default  = 3.5, floating point scalar
;       charsize  - character size, default = 2.0, floating point scalar
;       color     - color that the arrows and NE letters should be.  Default
;                    value is !P.COLOR
;       Data - if this keyword is set and nonzero, the input center (xcen,
;                 ycen) is understood to be in data coordinates
;       font - IDL vector font number (1-20) to use to display NE letters.
;                 For example, set font=13 to use complex italic font.
;       NotVertex - Normally (historically) the specified xcen,ycen indicated
;                   the position of the vertex of the figure.  If this
;                   keyword is set, the xcen,ycen coordinates refer to a sort
;                   of 'center of mass' of the figure.  This allows the
;                   figure to always appear with the area irregardless of
;                   the rotation angle.
;       Normal - if this keyword is set and nonzero, the input center 
;                (xcen,ycen) is taken to be in normalized coordinates.   The
;                default is device coordinates.
;       thick     - line thickness, default = 2.0, floating point scalar
; OUTPUTS:
;       none
; EXAMPLE:
;       Draw a weathervane at (400,100) on the currently active window, 
;       showing the orientation of the image associated with a FITS header, hdr
;
;       IDL> arrows, hdr, 400, 100
;
; METHOD:
;       Uses EXTAST to EXTract ASTrometry from the FITS header.   The 
;       directions of North and East are computed and the procedure
;       ONE_ARROW called to create the "weathervane".
;
; PROCEDURES USED:
;       EXTAST - Extract astrometry structure from FITS header
;       ONE_ARROW - Draw a labeled arrow	
;       ZPARCHECK
; REVISON HISTORY:
;       written by B. Boothman 2/5/86 
;       Recoded with new procedures ONE_ARROW, ONE_RAY.  R.S.Hill,HSTX,5/20/92
;       Added separate determination for N and E arrow to properly display
;         arrows irregardless of handedness or other peculiarities and added
;         /NotVertex keyword to improve positioning of figure. E.Deutsch 1/10/93
;       Added /DATA and /NORMAL keywords W. Landsman      July 1993
;       Recognize GSSS header    W. Landsman       June 1993
;       Added /FONT keyword W. Landsman           April 1995
;       Modified to work correctly for COLOR=0  J.Wm.Parker, HITC   1995 May 25
;       Work correctly for negative CDELT values   W. Landsman   Feb. 1996
;       Converted to IDL V5.0   W. Landsman   September 1997
;-

  On_error,2                            ;Return to caller

  if (N_params() LT 1) then begin 
    print,'Syntax - ' + $
             'ARROWS, hdr, [ xcen, ycen, ARROWLEN= , CHARSIZE=  COLOR= , /DATA
    print,'                        FONT=, /NORMAL, /NotVertex, THICK=  ]
    print,'         hdr - FITS header with astrometry
    return
  endif else zparcheck,'ARROWS',h,1,7,1,'FITS header array'

  if ( N_params() LT 3 ) then $
    read,'Enter x, y values for center of arrows: ',xcen,ycen

  if not keyword_set( THICK ) then    thick    = 2.0
  if not keyword_set( CHARSIZE ) then charsize = 2.0
  if not keyword_set( ARROWLEN ) then arrowlen = 3.5
  if (n_elements(COLOR) eq 0) then color = !P.COLOR
  if not keyword_set( NotVertex ) then NotVertex=0

;  Derive Position Angles for North and East separately
  extast, h, astr
  if strmid( astr.ctype[0], 5, 3)  EQ 'GSS' then begin
	h1 = h
	gsss_stdast, h1
	extast, h1, astr
  endif
  cd = astr.cd
  cdelt = astr.cdelt

  dRAdX = cdelt[0]*cd[0,0] & dRAdY = cdelt[1]*cd[0,1]
  dDECdX = cdelt[0]*cd[1,0] & dDECdY = cdelt[1]*cd[1,1]

  NPA = atan(dDECdY,dDECdX)*!radeg & if (NPA lt 0) then NPA=360+NPA
  EPA = atan(dRAdY,dRAdX)*!radeg & if (EPA lt 0) then EPA=360+EPA
  NPA = NPA-90 & EPA = EPA-90

;  Make arrows reasonable size depending on device

  arrowlen_dev = arrowlen*!D.y_ch_size
  arrowsize = [arrowlen_dev, arrowlen_dev/3.5, 35.0]  ; See one_arrow.pro

;  Adjust Center to 'Center of Mass' if NotVertex set
  if keyword_set( NORMAL) then begin
	newcen = convert_coord( xcen, ycen, /NORMAL, /TO_DEVICE)
        xcent = newcen[0]
        ycent = newcen[1]
  endif else if keyword_set( DATA) then begin
	newcen = convert_coord( xcen, ycen, /DATA, /TO_DEVICE)
        xcent = newcen[0]
        ycent = newcen[1]
  endif else begin
         xcent=xcen & ycent=ycen
  endelse 

 if NotVertex then begin
    RAnorm = sqrt( dRAdX^2 + dRAdY^2 )
    DECnorm = sqrt(dDECdX^2 + dDECdY^2 )
    xcent = xcen - (dRAdX+dDECdX)/2/RAnorm*arrowsize[0]
    ycent = ycen - (dRAdY+dDECdY)/2/DECnorm*arrowsize[0]
    endif

;  Draw arrows
  one_arrow, xcent, ycent,  90+NPA, 'N', font= font, $
    charsize=charsize, thick=thick, color=color, arrowsize=arrowsize
  one_arrow, xcent, ycent, 90+EPA, 'E', font = font, $
    charsize=charsize, thick=thick, color=color, arrowsize=arrowsize

  return
  end
