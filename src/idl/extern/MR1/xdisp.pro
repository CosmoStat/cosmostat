;
; Auto Save File For ./xdisp.pro
;
;  Tue Feb  7 17:39:25 MET 1995
;



; CODE MODIFICATIONS MADE ABOVE THIS COMMENT WILL BE LOST.
; DO NOT REMOVE THIS COMMENT: BEGIN HEADER


; $Id: xdisp.pro,v 1.2 1997/02/27 17:57:31 starck Exp $
;
;+
; NAME:
;	XDISP
;
; PURPOSE:
;	 XDISP is a widget program for image analysis. Several operations
;        can be done by using the mouse and press buttons:
;           QUIT: quit the application
;           LOAD: load an image. Image formats can be either fits or midas
;                 format
;           LUT: modify the LUT
;           PROFILE: examine lines or columns 
;           CURSOR: examine pixel values. If image format is fits and if
;                   the header contains astrometric position, then pixel
;                   position in the sky is given (Right ascension and
;                   declinaison.
;           HISTO: plot the histogram
;           INFO: print the min, max, mean, sigma of the image
;           FFT: computes the Fourier transform and display either the power 
;                spectrum, the phase, the real part or the imaginary part.
;           PAN: make a zoom. Zoom factor are 2,4,8,1/2,1/4,1/8
;
; CALLING SEQUENCE:
;	XDISP, Image, FitsHeader 
;
; INPUTS:
;	Image -- two dimensionnal array: image to analyse
;
; OPTIONAL INPUTS:
;	FitsHeader -- string array :	Fits header
;	
;
; EXAMPLE:
;               Image = readfits('data.fits', FitsHead)
;		XDISP, Image, FitsHead
;
;
; MODIFICATION HISTORY:
; 	Written by:	JL Starck 26/1/95 with widget builder
;       Add the contour and 3dim widget 9/2/95          
;
;-



; DO NOT REMOVE THIS COMMENT: END HEADER
; CODE MODIFICATIONS MADE BELOW THIS COMMENT WILL BE LOST.


; CODE MODIFICATIONS MADE ABOVE THIS COMMENT WILL BE LOST.
; DO NOT REMOVE THIS COMMENT: BEGIN PDMENU22



pro common_xdisp
COMMON  XDISP_BLOCK, WIDGET_SET, MAIN1, DRAW13, TEXT11,$
                 WDATA, Imag, FitsHeader, KFits, Zoom, ImagZoom, TabWin
end

PRO PDMENU22_Event, Event
COMMON  XDISP_BLOCK


  CASE Event.Value OF 


  'QUIT': BEGIN
     widget_control, MAIN1, /destroy
    END
  ENDCASE
END


; DO NOT REMOVE THIS COMMENT: END PDMENU22
; CODE MODIFICATIONS MADE BELOW THIS COMMENT WILL BE LOST.


; CODE MODIFICATIONS MADE ABOVE THIS COMMENT WILL BE LOST.
; DO NOT REMOVE THIS COMMENT: BEGIN PDMENU23




PRO PDMENU23_Event, Event
COMMON  XDISP_BLOCK

  CASE Event.Value OF 
  'LOAD': BEGIN
          case !version.os of
                  'vms': filt = '*.fits*'
                  'windows': filt = '*.fits'
;                  'MacOS': goto, notsup
                  else: filt = '*.fits'
               endcase

          FileName = pickfile(/READ, FILTER=filt, GET_PATH=GET_PATH)
          if strlen(FileName) GT 0 then BEGIN
            Format = 'unknown'
            if  !version.os EQ 'vms' then FileName = strlowcase(FileName)
            ; find the image format
            if rstrpos(FileName, '.d') GT 0  then Format = 'raw' else $
            if rstrpos(FileName, '.fit') GT 0 then Format = 'fits' else $
            if rstrpos(FileName, '.bdf') GT 0 then Format = 'midas' 
            KFits = 0
            CASE Format OF
              'midas': BEGIN
                     File = str_sep(FileName,'.bdf')
                     File = File[0]
                     MID_RD_IMAGE, File, Imag, NAXIS, NPIX
                       END
               'raw': Dat = rim(FileName)
               'fits': Imag = readfits(FileName, FitsHeader)
                else: BEGIN
                       print, 'xdisp: unable to read data, unknown format'
                       FitsHeader = 0
                       Format = 0
                       Imag = 0
                     END
            ENDCASE
            if (Format EQ 'midas') OR (Format EQ 'fits')  then BEGIN
               if type_of(FitsHeader) EQ 'STRING' then BEGIN
               extast,FitsHeader,astr,noparams
               if (noparams ge 0) then KFits = 1
               END
               Zoom = 1.
                erase
                vsize = size(Imag)
                Nl = vsize[2]
                Nc = vsize[1]
                N = Nl
                if N LT Nc then N = Nc
                while N*Zoom LT 512 do Zoom = Zoom * 2
                Nl = Nl*Zoom
                Nc = Nc*Zoom
                Imagzoom = congrid(Imag, Nc, Nl)
                tvscl, ImagZoom
                WIDGET_CONTROL, DRAW13, SET_DRAW_VIEW=[0,0]
             END
         endif
    END
  ENDCASE
END


; DO NOT REMOVE THIS COMMENT: END PDMENU23
; CODE MODIFICATIONS MADE BELOW THIS COMMENT WILL BE LOST.


; CODE MODIFICATIONS MADE ABOVE THIS COMMENT WILL BE LOST.
; DO NOT REMOVE THIS COMMENT: BEGIN PDMENU24


PRO PDMENU24_Event, Event


  CASE Event.Value OF 


  'LUT': BEGIN
    xloadct
    END
  ENDCASE
END


; DO NOT REMOVE THIS COMMENT: END PDMENU24
; CODE MODIFICATIONS MADE BELOW THIS COMMENT WILL BE LOST.


; CODE MODIFICATIONS MADE ABOVE THIS COMMENT WILL BE LOST.
; DO NOT REMOVE THIS COMMENT: BEGIN PDMENU36




PRO PDMENU36_Event, Event
COMMON  XDISP_BLOCK

vsize = size(imag)
if vsize[0] EQ 2 then BEGIN

  CASE Event.Value OF 


  'PRINT': BEGIN
     set_plot, 'PS' 
     device, filename='idl.ps'
     tvscl, ImagZoom
     device, /close
     set_plot, 'X'
    END
  ENDCASE
  ENDIF
END


; DO NOT REMOVE THIS COMMENT: END PDMENU36
; CODE MODIFICATIONS MADE BELOW THIS COMMENT WILL BE LOST.


; CODE MODIFICATIONS MADE ABOVE THIS COMMENT WILL BE LOST.
; DO NOT REMOVE THIS COMMENT: BEGIN PDMENU25




PRO PDMENU25_Event, Event
COMMON  XDISP_BLOCK
vsize = size(imag)
if vsize[0] EQ 2 then BEGIN


  CASE Event.Value OF 


  'GET.PROFILE': BEGIN
    str = 'Left mouse button to toggle between rows and columns. Left to Display and Right mouse button to Exit.'
      WIDGET_CONTROL,  TEXT11, SET_VALUE=str
    plot_profile, Imag, zoom=zoom, /order
    str = ' '
    WIDGET_CONTROL,  TEXT11, SET_VALUE=str
    END
  'GET.CURSOR': BEGIN
           Nx = (size(imag))[1]
           Ny = (size(imag))[2]
          str = [' Press left or center mouse button for new output line', $
                 '... right mouse button to exit.']
          WIDGET_CONTROL,  TEXT11, SET_VALUE=str
          if KFits EQ 0  then get_cursor, Imag, zoom=zoom $
          else  BEGIN
                get_cursor, FitsHeader, Imag, zoom=zoom
                END
          str = '  '
          WIDGET_CONTROL,  TEXT11, SET_VALUE=str
    END
  ENDCASE
ENDIF
END


; DO NOT REMOVE THIS COMMENT: END PDMENU25
; CODE MODIFICATIONS MADE BELOW THIS COMMENT WILL BE LOST.


; CODE MODIFICATIONS MADE ABOVE THIS COMMENT WILL BE LOST.
; DO NOT REMOVE THIS COMMENT: BEGIN PDMENU27


PRO PDMENU27_Event, Event
COMMON  XDISP_BLOCK

vsize = size(imag)
if vsize[0] EQ 2 then BEGIN
  CASE Event.Value OF 
  'PLOT.HISTO': BEGIN
              str = 'Right mouse button to Exit histogram.' 
          WIDGET_CONTROL,  TEXT11, SET_VALUE=str

           wsize = .75
           window,/free ,xs=wsize*640, ys=wsize*512,title='Histogram'  
           new_w = !d.window
           plothist, Imag
           print,'Right mouse button to Exit histogram.'
           Run = 1
	   while Run do begin
	       cursor,x,y,2,/dev	;Read position
               if !err eq 4 then  Run=0
           end
           wdelete, new_w
          str = ' '
          WIDGET_CONTROL,  TEXT11, SET_VALUE=str
    END
  'PLOT.CONTOURS': BEGIN
    a = !d.window
    xcont, Imag, window=9
    TabWin[9] = 1
    wset, TabWin[11]
    END
  'PLOT.3D VISU': BEGIN
    a = !d.window
    x3view, Imag
    wset, TabWin[11]
    END
  'PLOT.LOAD IN ANOTHER WINDOW.0': BEGIN
    window, 0, xsize=vsize[1]*Zoom, ysize=vsize[2]*Zoom
    tvscl, ImagZoom
    TabWin[0] = 1
    END
  'PLOT.LOAD IN ANOTHER WINDOW.1': BEGIN
    window, 1, xsize=vsize[1]*Zoom, ysize=vsize[2]*Zoom
    tvscl, ImagZoom
    TabWin[1] = 1
    END
  'PLOT.LOAD IN ANOTHER WINDOW.2': BEGIN
    window, 2, xsize=vsize[1]*Zoom, ysize=vsize[2]*Zoom
    tvscl, ImagZoom
    TabWin[2] = 1
    END
  'PLOT.LOAD IN ANOTHER WINDOW.3': BEGIN
    window, 3, xsize=vsize[1]*Zoom, ysize=vsize[2]*Zoom
    tvscl, ImagZoom  
    TabWin[3] = 1  
    END
  'PLOT.LOAD IN ANOTHER WINDOW.4': BEGIN
    window, 4, xsize=vsize[1]*Zoom, ysize=vsize[2]*Zoom
    tvscl, ImagZoom
    TabWin[4] = 1
    END
  'PLOT.LOAD IN ANOTHER WINDOW.5': BEGIN
    window, 5, xsize=vsize[1]*Zoom, ysize=vsize[2]*Zoom
    tvscl, ImagZoom
    TabWin[5] = 1
    END
  'PLOT.LOAD IN ANOTHER WINDOW.6': BEGIN
    window, 6, xsize=vsize[1]*Zoom, ysize=vsize[2]*Zoom
    tvscl, ImagZoom
    TabWin[6] = 1
    END
  'PLOT.LOAD IN ANOTHER WINDOW.7': BEGIN
    window, 7, xsize=vsize[1]*Zoom, ysize=vsize[2]*Zoom
    tvscl, ImagZoom
    TabWin[7] = 1
    END
  'PLOT.LOAD IN ANOTHER WINDOW.8': BEGIN
    window, 8, xsize=vsize[1]*Zoom, ysize=vsize[2]*Zoom
    tvscl, ImagZoom
    TabWin[8] = 1
    END
  'PLOT.LOAD IN ANOTHER WINDOW.9': BEGIN
    window, 9, xsize=vsize[1]*Zoom, ysize=vsize[2]*Zoom
    tvscl, ImagZoom
    TabWin[9] = 1
    END
  'PLOT.LOAD IN ANOTHER WINDOW.10': BEGIN
    window, 10, xsize=vsize[1]*Zoom, ysize=vsize[2]*Zoom
    tvscl, ImagZoom
    TabWin[10] = 1
    END
  ENDCASE
    wset, TabWin[11]
  ENDIF
END


; DO NOT REMOVE THIS COMMENT: END PDMENU27
; CODE MODIFICATIONS MADE BELOW THIS COMMENT WILL BE LOST.


; CODE MODIFICATIONS MADE ABOVE THIS COMMENT WILL BE LOST.
; DO NOT REMOVE THIS COMMENT: BEGIN PDMENU28




PRO PDMENU28_Event, Event
COMMON  XDISP_BLOCK

vsize = size(imag)
if vsize[0] EQ 2 then BEGIN


  CASE Event.Value OF 


  'INFO': BEGIN
    min = min(Imag)
    max = max(Imag)
    Mean = total(Imag) / N_ELEMENTS(Imag)
    Sigma = sqrt(total( (Imag-mean)^2) /N_ELEMENTS(Imag))
    str = 'min='+strcompress(string(min),/REMOVE_ALL)+ $
      ' max='+strcompress(string(max),/REMOVE_ALL)+ $
      ' mean='+strcompress(string(mean),/REMOVE_ALL)+ $
      ' sigma='+strcompress(string(sigma),/REMOVE_ALL)

    if Zoom NE 1 then str = str + ' zoom='+strcompress(string(zoom),/REMOVE_ALL)
      WIDGET_CONTROL,  TEXT11, SET_VALUE=str

    END
  ENDCASE
ENDIF
END


; DO NOT REMOVE THIS COMMENT: END PDMENU28
; CODE MODIFICATIONS MADE BELOW THIS COMMENT WILL BE LOST.


; CODE MODIFICATIONS MADE ABOVE THIS COMMENT WILL BE LOST.
; DO NOT REMOVE THIS COMMENT: BEGIN PDMENU7




PRO PDMENU7_Event, Event
COMMON  XDISP_BLOCK

vsize = size(Imag)
if vsize[0] EQ 2 then BEGIN

vsize = size(ImagZoom)
Nl = vsize[2]
Nc = vsize[1]

  CASE Event.Value OF 


  'FFT.POWER SPECTRUM': BEGIN
        Imag = power_spectrum(Imag)
    END
  'FFT.PHASE': BEGIN
    F = dft(Imag)
    Imag = atan(imaginary(f), float(f))
    END
  'FFT.REAL': BEGIN
    Imag = float(dft(Imag))
    END
  'FFT.IMAGINARY': BEGIN
    Imag = imaginary(dft(Imag))
    END
  ENDCASE
    erase
    if Zoom NE 1 then Imagzoom = congrid(Imag, Nc, Nl)  $
    else Imagzoom = Imag
    tvscl, Imagzoom
ENDIF
END


; DO NOT REMOVE THIS COMMENT: END PDMENU7
; CODE MODIFICATIONS MADE BELOW THIS COMMENT WILL BE LOST.


; CODE MODIFICATIONS MADE ABOVE THIS COMMENT WILL BE LOST.
; DO NOT REMOVE THIS COMMENT: BEGIN PDMENU9




PRO PDMENU9_Event, Event
COMMON  XDISP_BLOCK

vsize = size(imag)
if vsize[0] EQ 2 then BEGIN


  CASE Event.Value OF 


  'PAN.PAN 2': BEGIN
    fact = 2.
    END
  'PAN.PAN 4': BEGIN
   fact= 4.
    END
  'PAN.PAN 8': BEGIN
    fact = 8
    END
  'PAN.PAN 1/2': BEGIN
    fact = 0.5
    END
  'PAN.PAN 1/4': BEGIN
    fact = 0.25
    END
  'PAN.PAN 1/8': BEGIN
    fact = 0.125
    END
  ENDCASE
    vsize = size(ImagZoom)
    Nl = fix(vsize[2]*fact)
    Nc = fix(vsize[1]*fact)
    erase
    widget_control,/hourglass
    Imagzoom = congrid(Imag, Nc, Nl)
    zoom = zoom * fact
    tvscl, ImagZoom
END
END


; DO NOT REMOVE THIS COMMENT: END PDMENU9
; CODE MODIFICATIONS MADE BELOW THIS COMMENT WILL BE LOST.


; CODE MODIFICATIONS MADE ABOVE THIS COMMENT WILL BE LOST.
; DO NOT REMOVE THIS COMMENT: BEGIN MAIN1




PRO MAIN1_Event, Event
COMMON  XDISP_BLOCK

  WIDGET_CONTROL,Event.Id,GET_UVALUE=Ev

  CASE Ev OF 

  ; Event for PDMENU22
  'PDMENU22': PDMENU22_Event, Event
  ; Event for PDMENU23
  'PDMENU23': PDMENU23_Event, Event
  ; Event for PDMENU24
  'PDMENU24': PDMENU24_Event, Event
  ; Event for PDMENU36
  'PDMENU36': PDMENU36_Event, Event
  ; Event for PDMENU25
  'PDMENU25': PDMENU25_Event, Event
  ; Event for PDMENU27
  'PDMENU27': PDMENU27_Event, Event
  ; Event for PDMENU28
  'PDMENU28': PDMENU28_Event, Event
  ; Event for PDMENU7
  'PDMENU7': PDMENU7_Event, Event
  ; Event for PDMENU9
  'PDMENU9': PDMENU9_Event, Event
  'TEXT11': BEGIN
      Print, 'Event for TEXT11'
      END
  'DRAW13': BEGIN
            vsize = size(imag)
      if vsize[0] EQ 2 then BEGIN
           x = fix(event.x / zoom)
           y = fix(event.y / zoom)
           Nl = vsize[2]
           Nc = vsize[1]
           if (x ge 0 and x lt Nc) AND (y ge 0 and y lt Nl) then $
               Inten = Imag(x,y) else Inten = 0

           str =  '('+strcompress(string(x),/REMOVE_ALL)+','+$
                 strcompress(string(y),/REMOVE_ALL)+') = '+$
                 strcompress(string(Inten),/REMOVE_ALL)
           WIDGET_CONTROL,  TEXT11, SET_VALUE=str
          END
      END
  ENDCASE
END


; DO NOT REMOVE THIS COMMENT: END MAIN1
; CODE MODIFICATIONS MADE BELOW THIS COMMENT WILL BE LOST.



PRO xdisp, Dat, Hd, GROUP=Group 
COMMON  XDISP_BLOCK

  astrolib
  Zoom = 1
  Imag = 0
  TabWin = intarr(12)
  KFits  = 0

  if (size(Dat))[0] EQ 2 then Imag = Dat
  ImagZoom = Imag
  if N_PARAMS() EQ 2 then BEGIN
     if type_of(Hd) EQ 'STRING' then BEGIN
           extast,hd,astr,noparams
           if (noparams ge 0) then BEGIN
                KFits = 1
                FitsHeader = Hd
                END

          END
     END


  IF N_ELEMENTS(Group) EQ 0 THEN GROUP=0

  junk   = { CW_PDMENU_S, flags:0, name:'' }


  MAIN1 = WIDGET_BASE(GROUP_LEADER=Group, $
      /COLUMN, $
      MAP=1, $
      TITLE='XDISPLAY', $
      UVALUE='MAIN1', $
      XSIZE=600, $
      YSIZE=700)

  BASE3 = WIDGET_BASE(MAIN1, $
      ROW=2, $
      MAP=1, $
      UVALUE='BASE3')

  BASE6 = WIDGET_BASE(BASE3, $
      ROW=1, $
      MAP=1, $
      UVALUE='BASE6')

  MenuDesc1116 = [ $
      { CW_PDMENU_S,       2, 'QUIT' } $  ;      0

  ]


  PDMENU22 = CW_PDMENU( BASE6, MenuDesc1116, /RETURN_FULL_NAME, $
      UVALUE='PDMENU22')

  MenuDesc1117 = [ $
      { CW_PDMENU_S,       2, 'LOAD' } $  ;      0

  ]


  PDMENU23 = CW_PDMENU( BASE6, MenuDesc1117, /RETURN_FULL_NAME, $
      UVALUE='PDMENU23')

  MenuDesc1118 = [ $
      { CW_PDMENU_S,       2, 'LUT' } $  ;      0

  ]


  PDMENU24 = CW_PDMENU( BASE6, MenuDesc1118, /RETURN_FULL_NAME, $
      UVALUE='PDMENU24')

  MenuDesc1119 = [ $
      { CW_PDMENU_S,       2, 'PRINT' } $  ;      0

  ]


  PDMENU36 = CW_PDMENU( BASE6, MenuDesc1119, /RETURN_FULL_NAME, $
      UVALUE='PDMENU36')

  MenuDesc1120 = [ $
      { CW_PDMENU_S,       3, 'GET' }, $ ;        0
        { CW_PDMENU_S,       0, 'PROFILE' }, $ ;        1
        { CW_PDMENU_S,       2, 'CURSOR' } $  ;      2

  ]


  PDMENU25 = CW_PDMENU( BASE6, MenuDesc1120, /RETURN_FULL_NAME, $
      UVALUE='PDMENU25')

  MenuDesc1122 = [ $
      { CW_PDMENU_S,       3, 'PLOT' }, $ ;        0
        { CW_PDMENU_S,       0, 'HISTO' }, $ ;        1
        { CW_PDMENU_S,       0, 'CONTOURS' }, $ ;        2
        { CW_PDMENU_S,       0, '3D VISU' }, $ ;        3
        { CW_PDMENU_S,       3, 'LOAD IN ANOTHER WINDOW' }, $ ;        4
          { CW_PDMENU_S,       0, '0' }, $ ;        5
          { CW_PDMENU_S,       0, '1' }, $ ;        6
          { CW_PDMENU_S,       0, '2' }, $ ;        7
          { CW_PDMENU_S,       0, '3' }, $ ;        8
          { CW_PDMENU_S,       0, '4' }, $ ;        9
          { CW_PDMENU_S,       0, '5' }, $ ;       10
          { CW_PDMENU_S,       0, '6' }, $ ;       11
          { CW_PDMENU_S,       0, '7' }, $ ;       12
          { CW_PDMENU_S,       0, '8' }, $ ;       13
          { CW_PDMENU_S,       0, '9' }, $ ;       14
          { CW_PDMENU_S,       2, '10' } $  ;     15

  ]


  PDMENU27 = CW_PDMENU( BASE6, MenuDesc1122, /RETURN_FULL_NAME, $
      UVALUE='PDMENU27')

  MenuDesc1125 = [ $
      { CW_PDMENU_S,       2, 'INFO' } $  ;      0

  ]


  PDMENU28 = CW_PDMENU( BASE6, MenuDesc1125, /RETURN_FULL_NAME, $
      UVALUE='PDMENU28')

  MenuDesc1126 = [ $
      { CW_PDMENU_S,       3, 'FFT' }, $ ;        0
        { CW_PDMENU_S,       0, 'POWER SPECTRUM' }, $ ;        1
        { CW_PDMENU_S,       0, 'PHASE' }, $ ;        2
        { CW_PDMENU_S,       0, 'REAL' }, $ ;        3
        { CW_PDMENU_S,       2, 'IMAGINARY' } $  ;      4

  ]


  PDMENU7 = CW_PDMENU( BASE6, MenuDesc1126, /RETURN_FULL_NAME, $
      UVALUE='PDMENU7')

  MenuDesc1128 = [ $
      { CW_PDMENU_S,       3, 'PAN' }, $ ;        0
        { CW_PDMENU_S,       0, 'PAN 2' }, $ ;        1
        { CW_PDMENU_S,       0, 'PAN 4' }, $ ;        2
        { CW_PDMENU_S,       0, 'PAN 8' }, $ ;        3
        { CW_PDMENU_S,       0, 'PAN 1/2' }, $ ;        4
        { CW_PDMENU_S,       0, 'PAN 1/4' }, $ ;        5
        { CW_PDMENU_S,       2, 'PAN 1/8' } $  ;      6

  ]


  PDMENU9 = CW_PDMENU( BASE6, MenuDesc1128, /RETURN_FULL_NAME, $
      UVALUE='PDMENU9')


  TextVal1130 = [ '' ]
  TEXT11 = WIDGET_TEXT( BASE3,VALUE=TextVal1130, $
      UVALUE='TEXT11', $
      YSIZE=1)


  DRAW13 = WIDGET_DRAW( MAIN1, $
      BUTTON_EVENTS=1, $
      MOTION_EVENTS=1, $
      RETAIN=2, $
      UVALUE='DRAW13', $
      XSIZE=1024, $
      X_SCROLL_SIZE=512, $
      YSIZE=1024, $
      Y_SCROLL_SIZE=512, $
      /VIEWPORT_EVENTS)

  WIDGET_CONTROL, MAIN1, /REALIZE

  ; Get drawable window index

  COMMON DRAW13_Comm, DRAW13_Id
  WIDGET_CONTROL, DRAW13, GET_VALUE=DRAW13_Id

  vsize = size(Imag)
  if vsize[0] EQ 2 then BEGIN
             Nl = vsize[2]
             Nc = vsize[1]
             N = Nl
             if N LT Nc then N = Nc
             while (2*N)*Zoom LT 513 do Zoom = Zoom * 2
             Nl = Nl*Zoom
             Nc = Nc*Zoom
             Imagzoom = congrid(Imag, Nc, Nl)
             tvscl, ImagZoom
             END
  TabWin(11) = !d.window
  WIDGET_CONTROL, DRAW13, SET_DRAW_VIEW=[0,0]
  XMANAGER, 'MAIN1', MAIN1
END
