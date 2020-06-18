;+
; NAME:
;	X3D
;
; PURPOSE:
;	 X3D is a widget program for cube analysis. Several operations
;        can be done by using the mouse and press buttons:
;           QUIT: quit the application
;           Temporal Cut: An temporal cut is display
;           Horizontal Cut: An Horizontal cut is display
;           Vertical Cut: An Vertical Cut is display
;           Frame Number: When the user enter a value V and carriage return
;                         the image Cube(*,*,V) is displayed
;           Window Size: Define the number of elements plotted.
;           Slider Frame Number: When the user move the slider, 
;                        the corresponding frame is displayed
;
; CALLING SEQUENCE:
;	X3D, arg, [mask], GROUP=Group, noscale=noscale
;
; INPUTS:
;   arg      --   varies:             Cube (any 3-d numeric) OR
;                                     SCD  (with model)      OR
;                                     Structure (must have DATA & MASK fields)
;
; OPTIONAL INPUTS:
;   mask     --   byte array:         Only relevant if arg is a cube
; 	
; KEYED INPUTS:  
;   group    --   scalar:             widget group.
;   from     --   array :             BEGINning of displayed cubes.
;   to       --   array :             end of displayed cubes.
;   noscale  --   boolean:            if set don't scale undefined values
;
; EXAMPLE:
;  HELP,MY_CUBE 
;  MY_CUBE         INT       = Array(32, 32, 100)
;  X3D, my_cube,from=[10,20,50],to=[20,40,80] 
;
; ALGORITHM: created by widget builder.
;
; MODIFICATION HISTORY:
; 	Written by:	JL Starck 13/7/95 with widget builder
;       DL  15-May-1996 : Added keyed FROM and TO . 
;       S Guest, 21st  May 1996: xregistered, register as x3d;
;                                act only on button release;
;                                get_cbmsk allowing SCD & structure inputs.
;       S Guest, 17th July 1996: care with window selection.
;       DL  2x-Sep-1996 : Added LOAD LUT button . 
;                         Adding interval selection with mouse .
;       SG  11-Mar-1997 : undefined value scaling
;
;-

function scale_undef, im
image = im

type = type_code(image)
 CASE type OF
      1 : undef = 255b
      2 : undef = fix(-32768)
      3 : undef = -32768l
      4 : undef = -32768.
      5 : undef = -32768d
     endcase

igood = where (image ne undef, ngood)
ibad  = where (image eq undef, nbad)

if (nbad gt 0 and ngood gt 0) then begin
   mxgood = max(image[igood], min=mngood)

   ; set range. watch out for flat images
   if (mxgood ne mngood) then begin
       range = mxgood - mngood

   endif else begin
       case 1 of
           (mxgood gt 0): range = +mxgood
           (mxgood lt 0): range = -mxgood
            else        : range = 1
       endcase
   endelse

   ; replace undefined values with ones a bit smaller than everything else
   low = mngood - 0.05 * range
   image[ibad] = low
endif

return, image

end

;----------- common_x3d 
PRO common_x3d
common WIDGET_x3d, BASE_X3D, FIELDFrameNumber, DRAW_VECT_Id, DRAW_IMA_Id, $
                   SLIDER_FRAME, SLIDER_CUBE, FIELD_WINSIZE, $
                   fromto_flag, zoom_side_flag, bzoom_flag, zoom_state_flag, $
                   choose_plot_flag, AA, BB, linstyl_aa, linstyl_bb, $
                   mask_on, mask_id, Sx,Ex,Sy,Ey,St,Et
common GLOBAL_x3d, zoom, CubeX3D, MaskX3D, Ns, ImaZoom, Vector, Xc, Yc, $
                   TypeCut, WindowSize, Nx,Ny,Nz, SelectFrame, store_frame, ymax, cfrom, cto, $
                   SelectCube, scale
END

;-----------  display_image

PRO display_image
common WIDGET_x3d
common GLOBAL_x3d
;
 image = CubeX3d[*,*,SelectFrame]
;
 IF (mask_on) THEN BEGIN
    mask = MaskX3d[*,*,SelectFrame]
    i    = WHERE(mask GT 0, n)
    IF (n GT 0) THEN image[i] = la_undef(type_code(image))
 ENDIF

 if (scale) then image = scale_undef(image)
;
 ImaZoom = congrid ( reform(image), Ns, Ns)
 WSET,draw_ima_id
 TVSCL, ImaZoom
;
END 

;-----------  plot_cut

PRO plot_cut
 common WIDGET_x3d
 common GLOBAL_x3d
; Set 
 choose_plot_flag = 1 - bzoom_flag 
;
;;; PRINT,'choose_plot_flag:',choose_plot_flag
;;; PRINT,'bzoom_flag:',bzoom_flag
;;; PRINT,'zoom_state_flag:',zoom_state_flag
;
 CASE TypeCut OF
   'H': BEGIN
         IF ( (NOT bzoom_flag ) AND (NOT zoom_state_flag)) THEN BEGIN 
             Sx = max( [0, Xc-WindowSize/2] )
             Ex = min( [Nx-1, Xc+WindowSize/2] )
         ENDIF
          Vector = REFORM( cubeX3d[Sx:Ex, Yc, SelectFrame])
          W = INDGEN( (SIZE(vector))[1] ) + Sx
          txs = STRCOMPRESS(STRING(Sx), /REMOVE_ALL)
          txe = STRCOMPRESS(STRING(Ex), /REMOVE_ALL)
          ty  = STRCOMPRESS(STRING(Yc), /REMOVE_ALL)
          tt  = STRCOMPRESS(STRING(SelectFrame), /REMOVE_ALL)
          title='Horizontal cut: ('+txs+':'+txe+','+ty+','+tt+')'
          IF (mask_on) THEN Vmask = REFORM( maskX3d[Sx:Ex, Yc, SelectFrame])
         END
;
    'V': BEGIN
          IF ( (NOT bzoom_flag) AND (NOT zoom_state_flag)) THEN BEGIN 
             Sy = max( [0, Yc-WindowSize/2] )
             Ey = min( [Ny-1, Yc+WindowSize/2] )
          ENDIF
          Vector = REFORM( cubeX3d[Xc, Sy:Ey, SelectFrame])
          W = INDGEN( (SIZE(vector))[1] ) + Sy
          tt = STRCOMPRESS(STRING(SelectFrame), /REMOVE_ALL)
          tx = STRCOMPRESS(STRING(Xc), /REMOVE_ALL)
          tys = STRCOMPRESS(STRING(Sy), /REMOVE_ALL)
          tye = STRCOMPRESS(STRING(Ey), /REMOVE_ALL)
          title='Vertical cut: (' +tx + ',' + tys + ':' + tye + ',' + tt + ')'
          IF (mask_on) THEN Vmask = reform( maskX3d[Xc, Sy:Ey, SelectFrame])
         END
;
    'T': BEGIN
         IF ( (NOT bzoom_flag) AND (NOT zoom_state_flag) ) THEN BEGIN 
           IF ( fromto_flag ) THEN BEGIN
              St = cFrom[SelectCube-1]
              Et = cTo[SelectCube-1]
           ENDIF ELSE BEGIN
              St = MAX( [0, SelectFrame-WindowSize/2] )
              Et = MIN( [Nz-1, SelectFrame+WindowSize/2] )
           ENDELSE
         ENDIF
;
          Vector = REFORM( cubeX3d[Xc, Yc, St:Et])
          W = INDGEN( (size(vector))[1] ) + St
          tx = STRCOMPRESS(STRING(Xc), /REMOVE_ALL)
          ty = STRCOMPRESS(STRING(Yc), /REMOVE_ALL)
          tts = STRCOMPRESS(STRING(St), /REMOVE_ALL)
          tte = STRCOMPRESS(STRING(Et), /REMOVE_ALL)
          title='Temporal cut: ('+tx +',' + ty +','+tts+':'+tte+')'
;
          IF (mask_on) THEN Vmask = reform( maskX3d[Xc, Yc, St:Et])
         END
   ENDCASE
;
 WSET, DRAW_VECT_Id
 p = !p.multi
 !p.multi = 0
 PLOT, w, Vector, TITLE=title, YSTYLE=18
 !p.multi = p
;
 IF (mask_on) THEN BEGIN
    i = WHERE(Vmask NE 0, n)
    IF (n gt 0) THEN OPLOT, w[i], Vector[i], psym=1
 ENDIF
;
;;; PRINT,'>>>> plot_cut :  St :',St
; Get last Y range for plot
 YY = [!Y.crange[0],!Y.crange[1]]
 ymax = MAX(YY)
END

;-----------------------------
PRO X3D_Event, Event
common WIDGET_x3d
common GLOBAL_x3d

  WIDGET_CONTROL,Event.Id,GET_UVALUE=Ev

  CASE Ev OF 

  'BGROUP_CUT': BEGIN
       CASE Event.Value OF
       0: BEGIN
           TypeCut = 'T'
          END
       1: BEGIN
           TypeCut = 'H'
          END
       2: BEGIN
           TypeCut = 'V'
          END
       ELSE: Message, 'Unknown button pressed'
       ENDCASE
       plot_cut
      END
;
  'FIELDFrameNumber': BEGIN
;
       WIDGET_CONTROL, FIELDFrameNumber, GET_VALUE=SelectFrame
;
       IF (fromto_flag EQ 0) THEN BEGIN 
          IF (SelectFrame GT Nz-1) THEN BEGIN
               SelectFrame =  Nz-1
               WIDGET_CONTROL, FIELDFrameNumber, SET_VALUE=SelectFrame
;;;               PRINT,'Nz :',Nz
          ENDIF
;
       ENDIF ELSE BEGIN
;;;          PRINT,'fromto_flag :',fromto_flag
          IF ( NOT bzoom_flag ) THEN BEGIN
              IF SelectFrame LT cFROM(SelectCube-1) THEN SelectFrame = cFROM(SelectCube-1)      
              IF SelectFrame GT cTO(SelectCube-1) THEN SelectFrame = cTO(SelectCube-1)      
          ENDIF ELSE BEGIN
              set_plot_range
          ENDELSE
          WIDGET_CONTROL, FIELDFrameNumber, SET_VALUE=SelectFrame
       ENDELSE
;
;;;       PRINT,'SelectFrame :',SelectFrame
       WIDGET_CONTROL,  SLIDER_FRAME, SET_VALUE=SelectFrame
;
       display_image
       plot_cut
      END
;
  'FIELD_WINSIZE': BEGIN
       WIDGET_CONTROL, FIELD_WINSIZE, GET_VALUE=chWindowSize
       WindowSize=FIX(chWindowSize)
       IF (fromto_flag) THEN BEGIN
;         help, WindowSize
;
;         HELP,cTO 
;;;         PRINT,'>> size cTO   :',SIZE(cTO)
;;;         PRINT,'>> SelectCube :',SelectCube
;;;         PRINT,'>> WindowSize :',WindowSize
          DIFF = cTO[SelectCube-1]-cFROM[SelectCube-1]+1
;;;         PRINT,'MAX :',MAX(WindowSize,DIFF) 
;
          WindowSize=MAX(WindowSize,DIFF) 
;;;         PRINT,'>> WindowSize :',WindowSize
       ENDIF
;
       plot_cut
      END
;
  'SLIDER_FRAME': BEGIN
       WIDGET_CONTROL,  SLIDER_FRAME, GET_VALUE=SelectFrame
       WIDGET_CONTROL, FIELDFrameNumber, SET_VALUE=SelectFrame
       display_image
       IF (NOT zoom_state_flag) THEN plot_cut
       choose_plot_flag = 0
      END
;
  'SLIDER_CUBE' : BEGIN
       FromTo_flag = 1
       WIDGET_CONTROL,  SLIDER_CUBE, GET_VALUE=SelectCube   
       SelectFrame = cFROM[SelectCube-1]
       WIDGET_CONTROL, FIELD_WINSIZE, GET_VALUE=chWindowSize
       WindowSize=FIX(chWindowSize)
       ImaZoom = CONGRID( REFORM(CubeX3d[*,*,SelectFrame] ), Ns, Ns)
       WSET,DRAW_IMA_Id
;
       St = cFROM[SelectCube-1]
       Et = cTO[SelectCube-1]
;
       WIDGET_CONTROL, FIELDFrameNumber, SET_VALUE=SelectFrame                              
       WIDGET_CONTROL, SLIDER_FRAME, SET_SLIDER_MIN=Et  ;;; I don't know why, but it works only like this (LD) !!!!
       WIDGET_CONTROL, SLIDER_FRAME, SET_SLIDER_MIN=St
       WIDGET_CONTROL, SLIDER_FRAME, SET_SLIDER_MAX=Et
       WIDGET_CONTROL, SLIDER_FRAME, SET_VALUE=SelectFrame
       display_image
       plot_cut
;
      END
;
  'DRAW_VECT': BEGIN

;       Zoom should not be displayed after this time .
        zoom_state_flag = 0

;       Zoom range selection :
        IF (Event.press) THEN BEGIN
              CC = CONVERT_COORD(Event.X,Event.Y,/DEVICE,/TO_DATA)
              WSET, DRAW_VECT_Id
              IF ( choose_plot_flag ) THEN BEGIN
                 IF ( zoom_side_flag ) THEN BEGIN
                    BB = CC
                    zoom_side_flag = 1 - zoom_side_flag
                    linstyl_aa = 2
                    linstyl_bb = 0
                 ENDIF ELSE BEGIN
                    AA = CC
                    zoom_side_flag = 1 - zoom_side_flag
                    linstyl_aa = 0
                    linstyl_bb = 2
                 ENDELSE
              ENDIF
;;
              plot_cut
;;
              !p.linestyle = linstyl_aa 
              OPLOT,[AA[0],AA[0]],[-ymax,ymax]
              !p.linestyle = linstyl_bb 
              OPLOT,[BB[0],BB[0]],[-ymax,ymax]
              !p.linestyle = 0
              AB=[AA[0],BB[0]]
              MN=MIN(AB)
              MX=MAX(AB)
              XYOUTS,150,10,'ymax '+strn(ymax),/DEVICE                   
;;
              IF (NOT zoom_side_flag) THEN BEGIN
                 XYOUTS,0,10,'ZOOM ',/DEVICE                   
                 XYOUTS,30,10,'FROM : '+strn(MN[0])+'<<<',/DEVICE    
                 XYOUTS,30,0,'TO : '+strn(MX[0]),/DEVICE       
              ENDIF ELSE BEGIN
                 XYOUTS,0,10,'ZOOM ',/DEVICE                   
                 XYOUTS,30,10,'FROM : '+strn(MN(0)),/DEVICE    
                 XYOUTS,30,0,'TO : '+strn(MX[0])+'<<<',/DEVICE       
              ENDELSE
;;
              WIDGET_CONTROL, SLIDER_FRAME, SET_SLIDER_MIN=Et  ;;; I don't know why, but it works only like this !!!!
              WIDGET_CONTROL, SLIDER_FRAME, SET_SLIDER_MIN=St
              WIDGET_CONTROL, SLIDER_FRAME, SET_SLIDER_MAX=Et
              WIDGET_CONTROL, SLIDER_FRAME, SET_VALUE=St 
;;
         ENDIF
;
       END
;
  'DRAW_IMA': BEGIN
      IF (event.type eq 1) THEN BEGIN   ; act on release only
           x = FIX(event.x / zoom)
           IF !order EQ 1 THEN y = Ny - fix(event.y / zoom) - 1 $
           ELSE y = FIX(event.y / zoom)
           Xc = x
           Yc = y
           tx = STRCOMPRESS(string(Xc), /REMOVE_ALL)
           ty = STRCOMPRESS(string(Yc), /REMOVE_ALL)
           tt = STRCOMPRESS(string(SelectFrame), /REMOVE_ALL)
           PRINT,'Cube('+tx +',' + ty +','+tt+') = ', cubeX3d[Xc,Yc,SelectFrame]
;
           WSET,DRAW_IMA_ID
           display_image
           XYOUTS,event.x,event.y,'+',/device
           plot_cut
;
        ENDIF
        END
;
  'BUTTON_NEXT': BEGIN
         SelectFrame = SelectFrame + 1                                                        
;                                                                                             
         IF ( SelectFrame GT Et ) THEN SelectFrame = Et
;                                                                                             
         WIDGET_CONTROL, FIELDFrameNumber, SET_VALUE=SelectFrame                              
         WIDGET_CONTROL,  SLIDER_FRAME, SET_VALUE=SelectFrame                                 
         display_image                                                                        
         IF (NOT zoom_state_flag) THEN plot_cut
      END
;
  'BUTTON_PREVIOUS': BEGIN
          SelectFrame = SelectFrame - 1                                                    
         IF ( SelectFrame LT St ) THEN SelectFrame = St
;;
          WIDGET_CONTROL, FIELDFrameNumber, SET_VALUE=SelectFrame                          
          WIDGET_CONTROL,  SLIDER_FRAME, SET_VALUE=SelectFrame                             
          display_image                                                                    
          IF (NOT zoom_state_flag) THEN plot_cut
      END
;
  'BUTTON_ALL': BEGIN
                    FromTo_flag=0                                                  
                    St = 0                                                         
                    Et = Nz-1                                                      
;;;                    PRINT,"BUTTON_ALL: Et",Et                                      
;;;                    PRINT,"BUTTON_ALL: St",St                                      
                    SelectFrame=1                                                  
                    plot_cut                                                       
                    WIDGET_CONTROL, SLIDER_FRAME, SET_SLIDER_MIN=St                
                    WIDGET_CONTROL, SLIDER_FRAME, SET_SLIDER_MAX=Et                
;;;                    PRINT,'cFROM(',SelectCube-1,') :',cFROM[SelectCube-1]          
                    WIDGET_CONTROL,  SLIDER_FRAME, SET_VALUE=cFROM[SelectCube-1]   
                    SelectFrame=cFROM[SelectCube-1]                                                   
                    display_image                                                  
                END                                                               

  'BUTTON_MASK': BEGIN
      mask_on = Event.Select
      display_image
      plot_cut
   END

  'BUTTON_LUT': BEGIN
         XLOADCT
      END

  'BUTTON_ZOOM': BEGIN
       store_frame = selectFrame
       IF ( (AA[0] NE 0) OR (BB[0] NE 0) ) THEN BEGIN 
             set_plot_range
             bzoom_flag=1
             plot_cut
             bzoom_flag=0
             zoom_state_flag = 1
             WIDGET_CONTROL, FIELDFrameNumber, SET_VALUE=St
             WIDGET_CONTROL, SLIDER_FRAME, SET_SLIDER_MIN=Et  ;;; I don't know why, but it works only like this !!!!
             WIDGET_CONTROL, SLIDER_FRAME, SET_SLIDER_MIN=St
             WIDGET_CONTROL, SLIDER_FRAME, SET_SLIDER_MAX=Et
             WIDGET_CONTROL, SLIDER_FRAME, SET_VALUE=St
;;;;;;;;
             display_image
;
             WSET, DRAW_VECT_Id
             XYOUTS,0,10,'CLICK ON THIS FRAME TO GET BACK PREVIOUS PLOT  .',/DEVICE                   
;
       ENDIF ELSE BEGIN
             bzoom_flag=0
             zoom_state_flag=0
             WSET, DRAW_VECT_Id
             XYOUTS,0,10,'PLEASE SELECT ZOOM RANGE ON THIS PLOT. ',/DEVICE                   
       ENDELSE
;
      END

  'BUTTON_QUIT': BEGIN
         WIDGET_CONTROL,/DESTROY, BASE_X3D 
      END
;
  ENDCASE
;
END

;-------------------------------------------------
PRO set_plot_range
common WIDGET_x3d
common GLOBAL_x3d
;
    CASE TypeCut OF                                                                            
       'H': BEGIN                                                                              
             Sx = AA[0]                                                                        
             Ex = BB[0]                                                                        
;;;             PRINT,'Sx : ',Sx                                                                  
;;;             PRINT,'Ex : ',Ex                                                                  
             dd = Ex-Sx                                                                        
             IF (dd LT 0) THEN BEGIN                                                           
                dd = Ex                                                                        
                Ex = Sx                                                                        
                Sx = dd                                                                        
             ENDIF                                                                             
           END                                                                                 
       'V': BEGIN                                                                              
              Sy = AA[0]                                                                       
              Ey = BB[0]                                                                       
;;;              PRINT,'Sy : ',Sy                                                                 
;;;              PRINT,'Ey : ',Ey                                                                 
              dd = Ey-Sy                                                                       
              IF (dd LT 0) THEN BEGIN                                                          
                 dd = Ey                                                                       
                 Ey = Sy                                                                       
                 Sy = dd                                                                       
              ENDIF                                                                            
            END                                                                                
       'T': BEGIN                                                                              
              St = AA[0]                                                                       
              Et = BB[0]                                                                       
;;;              PRINT,'St : ',St                                                                 
;;;              PRINT,'Et : ',Et                                                                 
              dd=Et-St                                                                         
              IF (dd LT 0) THEN  BEGIN                                                         
                 dd = Et                                                                       
                 Et = St                                                                       
                 St = dd                                                                       
              ENDIF                                                                            
            END                                                                                
    ENDCASE                                                                                    
;
    IF (SelectFrame LT St) THEN SelectFrame = St
    IF (SelectFrame GT Et) THEN SelectFrame = Et
;
END
;-------------------------------------------------

PRO X3D, arg1, arg2, GROUP=Group, FROM=from, TO=to, noscale=noscale
common WIDGET_x3d
common GLOBAL_x3d
;
; has common blocks, so only one copy allowed
  IF xregistered('x3d') THEN return

 IF N_PARAMS() LT 1 THEN BEGIN
   PRINT, 'CALLING SEQUENCE: ', $ 
    'X3D, arg1, arg2, GROUP=Group, FROM=from, TO=to, noscale=noscale'
    GOTO, CLOSING
 ENDIF

; get cube & mask, whatever the input was
;  get_cbmsk, arg1, arg2, cubeX3d, maskX3d, ack=ok
;  IF (NOT ok) THEN return
cubeX3d = arg1
IF N_PARAMS() LT 2 then begin
    maskX3d = arg1
    maskX3d(*) = 0
    END
    
  scale = not keyword_set(noscale)

; extract dimensions
  szCube  = size(cubeX3d)
  Nx      = szCube[1]
  Ny      = szCube[2]
  Nz      = szCube[3]

; the mask may not be defined
  mask_on = keyword_set(maskX3d)
;
  Nfrom = (size(from))
  Nto   = (size(to))
;;;  PRINT,'>>> Nfrom :',Nfrom
;;;  PRINT,'>>> Nto   :',Nto
;
; Parameters for zoom on plots :
  Sx=0
  Ex=0
  Sy=0
  Ey=0
  St=0
  Et=0
;
  zoom_side_flag = 0    ; flag for beginning or end of zoom interval definition
  bzoom_flag = 0        ; flag for zoom_button action
  zoom_state_flag = 0   ; flag for zoom status 
  store_frame = 0
  AA = INTARR(3)
  BB = INTARR(3)
  linstyl_aa = 2
  linstyl_bb = 0
; Get last Y range for plot
  YY = [!Y.RANGE(0),!Y.RANGE(1)]
  ymax = MAX(YY)
;
  xsize = 600 + 10 * FIX (ALOG10( NZ )) 
  IF ( keyword_set(from) and keyword_set(to) ) THEN BEGIN
    IF (Nfrom[0] NE Nto[0]) THEN BEGIN
       PRINT,"ARRAYS From[] and To[] should have same size ."
       return
    ENDIF 
    cfrom = from
    cto = to
    selectCube = 1
    fromto_flag=1
    SelectFrame = cfrom[0]
    ysize = 805
  ENDIF ELSE BEGIN
    cfrom=INTARR(1)
    cto=INTARR(1)
    fromto_flag=0
    SelectFrame = Nz / 2
    selectCube = 1
    ysize = 800
  ENDELSE

;
;;;  PRINT,">>> SelectFrame :",SelectFrame
  nn = Nx
  IF nn lt Ny THEN nn = Ny
  zoom = 320. / nn
  Xc = Nx / 2
  Yc = Ny / 2
  Ns = 320
  ImaZoom = CONGRID( REFORM(CubeX3d[*,*,SelectFrame] ), Ns, Ns)
  MaxWindowSize = 2*MAX([Nz, 32])
  WindowSize = MaxWindowSize

  TypeCut = 'T'

  IF N_ELEMENTS(Group) EQ 0 THEN GROUP=0

  junk   = { CW_PDMENU_S, flags:0, name:'' }

  BASE_X3D = WIDGET_BASE(GROUP_LEADER=Group, $
      COLUMN=1, $
      MAP=1, $
      TITLE='Cube Analysis', $
      UVALUE='BASE_X3D', $
      XSIZE=xsize, $
      YSIZE=ysize)

  BASE_1 = WIDGET_BASE(BASE_X3D, $
      ROW=1, $
      MAP=1, $
      SPACE=30, $
      UVALUE='BASE_1', $
      /FRAME)

  BASE_1_2 = WIDGET_BASE(BASE_1, $
      COLUMN=1, $
      MAP=1, $
      XPAD=10, $
      YPAD=25, $
      UVALUE='BASE_1_2' ,$
      /FRAME)

  BUTTON_LUT = WIDGET_BUTTON( BASE_1_2, $
      UVALUE='BUTTON_LUT', $
      VALUE='LOAD LUT')

  BUTTON_QUIT = WIDGET_BUTTON( BASE_1_2, $
      UVALUE='BUTTON_QUIT', $
      VALUE='Quit')

  xbase   = WIDGET_BASE(base_1_2, /exclusive)
;
  mask_id = WIDGET_BUTTON(xbase, uvalue='BUTTON_MASK', value='Mask')
  IF (not mask_on) THEN WIDGET_CONTROL, mask_id, SENSITIVE=0
  mask_on = 0             ; off by default

  BASE_1_3 = WIDGET_BASE(BASE_X3D, $
       ROW=1, $
       MAP=1, $                                         
       UVALUE='BASE_1_3' ,$
       /FRAME)

  IF ( ( Nfrom[0] EQ 1 ) AND ( Nto[0] EQ 1) ) THEN BEGIN 

;;;     PRINT,'Nfrom[1]:',Nfrom[1]
     SLIDER_CUBE = WIDGET_SLIDER( BASE_1_3, $ ;   BASE_1_3
        MINIMUM=1, $
        MAXIMUM=Nfrom(1), $
        XSIZE=320, $
        XOFFSET=0, $
        TITLE='SubCube Number', $
        UVALUE='SLIDER_CUBE', $
        VALUE=1 )
;
     BUTTON_ALL = WIDGET_BUTTON( BASE_1_3, $
        UVALUE='BUTTON_ALL', $
        VALUE='ALL FRAMES')

  ENDIF
;
  BUTTON_ZOOM = WIDGET_BUTTON( BASE_1_3, $
        UVALUE='BUTTON_ZOOM', $
        VALUE='ZOOM')

  ValBtns_Cut = [ $
    'Temporal cut', $
    'Horizontal cut', $
    'Vertical cut' ]

  BGROUP_CUT = CW_BGROUP( BASE_1, ValBtns_Cut, $
      ROW=3, $
      EXCLUSIVE=1, $
      XPAD=50, $
      SET_VALUE=0, $
      UVALUE='BGROUP_CUT')

  BASE_1_1 = WIDGET_BASE(BASE_1, $
      ROW=2, $
      MAP=1, $
      UVALUE='BASE_1_1', $
      /FRAME)

  FieldValFrame = [ STRING(SelectFrame) ]
;  HELP,FieldValFrame

  FIELDFrameNumber = CW_FIELD(BASE_1_1,VALUE=FieldValFrame, $
      ROW=1, $
      INTEGER=1, $
      TITLE='Frame Number', $
      UVALUE='FIELDFrameNumber', $
      XSIZE=4, /RETURN_EVENTS)

  FieldValWindSize = [ string(WindowSize) ]
  FIELD_WINSIZE = CW_FIELD( BASE_1_1,VALUE=FieldValWindSize, $
      ROW=1, $
      INTEGER=1, $
      TITLE='Window Size', $
      UVALUE='FIELD_WINSIZE', $
      XSIZE=4, /RETURN_EVENTS)

  BASE_1_3 = WIDGET_BASE(BASE_1, $
      COLUMN=1, $
      MAP=1, $
      UVALUE='BASE_1_3')
;

  DRAW_VECT = WIDGET_DRAW( BASE_X3D, $
      BUTTON_EVENTS=1, $
      FRAME=4, $
      RETAIN=2, $
      UVALUE='DRAW_VECT', $
      XSIZE=500, $
      YSIZE=150)

  BASE_SHOW= WIDGET_BASE(BASE_X3D, $  
                 ROW=1, $ 
                 YSIZE=400)

  BASE_BUTTONS = WIDGET_BASE(BASE_SHOW, $
      COLUMN=1, $
      XPAD=5, $
      MAP=1, $
      UVALUE='BASE_BUTTONS',$
      /FRAME )

  BUTTON_NEXT = WIDGET_BUTTON( BASE_BUTTONS, $
      UVALUE='BUTTON_NEXT', $
      VALUE='Next frame')

  BUTTON_PREVIOUS = WIDGET_BUTTON( BASE_BUTTONS, $
      UVALUE='BUTTON_PREVIOUS', $
      VALUE='Previous frame')

  SLIDER_FRAME = WIDGET_SLIDER( BASE_SHOW, $ ;BASE_1_2
      MAXIMUM=Nz-1, $
      MINIMUM=0, $
      YSIZE=320, $
      XOFFSET=0, $
      TITLE=' ', $
      UVALUE='SLIDER_FRAME', $
      VALUE=SelectFrame, $ 
      /VERTICAL )

  BASE_IMA = WIDGET_BASE(BASE_SHOW, $  
      ROW=1, $
      XPAD=5, $
      MAP=1, $
      UVALUE='BASE_IMA', $
      YSIZE=600)

  DRAW_IMA = WIDGET_DRAW( BASE_IMA, $
      BUTTON_EVENTS=1, $
      RETAIN=2, $
      UVALUE='DRAW_IMA', $
      X_SCROLL_SIZE=320, $
      Y_SCROLL_SIZE=320, $
      XSIZE=320, $
      YSIZE=320)
  ;
  WIDGET_CONTROL, BASE_X3D, /REALIZE

  ; Get drawable window index
  WIDGET_CONTROL, DRAW_VECT, GET_VALUE=DRAW_VECT_Id 

  ; Get drawable window index
  WIDGET_CONTROL, DRAW_IMA, GET_VALUE=DRAW_IMA_Id

  display_image
  plot_cut

  XMANAGER, 'x3d', BASE_X3D
;  CubeX3d =0

CLOSING:
;
END
