;+
; NAME:
;     XMOVIE
;
; PURPOSE:
;     This program is a simplified version of XINTERANIMATE. It is written
;     to illustrate the proper way to write an animation loop in a widget
;     program using the WIDGET_TIMER functionality.
;
; CATEGORY:
;     Widgets.
;
; CALLING SEQUENCE:
;
;      XMOVIE, image3d
;
; INPUTS:
;     image3d: A three-dimensional image array. The animation occurs over
;              the third dimension.
;
; KEYWORD PARAMETERS:
;     GROUP:   The group leader of the program. When the group leader dies,
;              this program dies as well.
;
;     TITLE:   The window title of the program. The default is "Animation
;              Example...".
;
; COMMON BLOCKS:
;     None.
;
; SIDE EFFECTS:
;     None.
;
; EXAMPLE:
;     To open the abnormal heart data and animate it, type:
;
;        filename = FILEPATH(SUBDIR=['examples', 'data'], 'abnorm.dat')
;        OPENR, lun, filename, /GET_LUN
;        data = BYTARR(64, 64, 15)
;        READU, lun, data
;        FREE_LUN, lun
;        data = REBIN(data, 256, 256, 15)
;
;        XMOVIE, data
;
; MODIFICATION HISTORY:
;      Written by: David Fanning, June 96.
;      Added slider for controlling animation speed. 30 June 99.
;-


PRO XMOVIE_DELAY_EVENTS, event

   ; Get the INFO structure.

WIDGET_CONTROL, event.top, GET_UVALUE=info, /NO_COPY

   ; Set the delay.

info.delay = (100 - event.value) / 200.0 ; Max delay 0.5 seconds.

      ; Put INFO structure back.

WIDGET_CONTROL, event.top, SET_UVALUE=info, /NO_COPY


END ; ***************************** of XMOVIE_DELAY_EVENTS ***********



PRO XMOVIE_EVENT, event

   ; Get the INFO structure.

WIDGET_CONTROL, event.top, GET_UVALUE=info, /NO_COPY

   ; What kind of event is this?

eventName = TAG_NAMES(event, /STRUCTURE_NAME)

   ; Code for BUTTON events.

IF eventName EQ 'WIDGET_BUTTON' THEN BEGIN

   CASE event.id OF

      info.start: BEGIN

            ; If stop flag is 0, then display next frame.

         IF info.stopflag EQ 0 THEN BEGIN
            WSET, info.wid
            TVSCL, info.data(*,*,info.thisFrame)

               ; Write the frame number on the image.

            XYOUTS, 0.05, 0.05, /NORMAL, 'Frame ' + STRTRIM(info.thisFrame,2)
         ENDIF

            ; Update frame number.

         info.thisFrame = info.thisFrame + 1
         IF info.thisFrame GT (info.nframes - 1) THEN info.thisFrame = 0

            ; Set a timer event to get back into this event handler.

         WIDGET_CONTROL, event.ID, TIMER=0

            ; Update stop flag.

         info.stopflag = 0
      END ; of START button CASE.

      info.stop: BEGIN
         info.stopflag = 1
         info.thisFrame = info.thisFrame - 1  ; Subtract a frame.
         END

      info.quit: WIDGET_CONTROL, event.top, /DESTROY
   ENDCASE

      ; Put INFO structure back.

   IF WIDGET_INFO(event.top, /VALID_ID) THEN $
      WIDGET_CONTROL, event.top, SET_UVALUE=info, /NO_COPY
   RETURN
ENDIF ; of BUTTON code

   ; Code for TIMER events.

IF eventName EQ 'WIDGET_TIMER' THEN BEGIN

   IF info.stopflag EQ 0 THEN BEGIN
      WSET, info.wid

         ; Delay for a moment.

      Wait, info.delay

         ; Display the image.

      TVSCL, info.data(*,*,info.thisFrame)

         ; Write the frame number on the image.

      XYOUTS, 0.05, 0.05, /NORMAL, 'Frame ' + STRTRIM(info.thisFrame,2)

         ; Update frame number.

      info.thisFrame = info.thisFrame + 1
      IF info.thisFrame GT (info.nframes - 1) THEN info.thisFrame = 0

         ; Set a timer event to get back into this event handler.

      IF info.stopflag EQ 0 THEN WIDGET_CONTROL, event.ID, TIMER=0
   ENDIF

      ; Put INFO structure back.

   WIDGET_CONTROL, event.top, SET_UVALUE=info, /NO_COPY
ENDIF ; of TIMER code

END ; ***************************** of XMOVIE_EVENT ***********************



PRO XMOVIE, data, GROUP=group, TITLE=title

ON_ERROR, 1

   ; Read some 3D example data for animation if none passed in.

IF N_PARAMS() EQ 0 THEN BEGIN
   filename = FILEPATH(SUBDIR=['examples', 'data'], 'abnorm.dat')
   OPENR, lun, filename, /GET_LUN
   data = BYTARR(64, 64, 15)
   READU, lun, data
   FREE_LUN, lun
   data = REBIN(data, 256, 256, 15)
ENDIF

IF N_ELEMENTS(title) EQ 0 THEN title='Animation Example...'

  ; Collect information about the dimensions of the data.

dataSize = SIZE(data)
IF dataSize(0) NE 3 THEN MESSAGE, 'Data parameter must be 3D.'

   ; Create the widgets.

tlb = WIDGET_BASE(COLUMN=1, TITLE=title)
start = WIDGET_BUTTON(tlb, Value='Start Animation')
stop = WIDGET_BUTTON(tlb, Value='Stop Animation')
drawID = WIDGET_DRAW(tlb, XSIZE=dataSize(1), YSIZE=dataSize(2))
delayID = WIDGET_SLIDER(tlb, Max=100, Min=0, /Suppress_Value, $
   Title='Frame Delay', Value=50, Event_Pro='XMovie_Delay_Events')
quit = WIDGET_BUTTON(tlb, Value='Exit Program')

   ; Realize and make draw widget active graphics window.

WIDGET_CONTROL, tlb, /REALIZE
WIDGET_CONTROL, drawID, GET_VALUE=wid
WSET, wid

   ; Load color table and first frame.

Device, Decomposed=0
LOADCT, 3
TVSCL, data(*,*,0)
XYOUTS, 0.05, 0.05, /NORMAL, 'Frame 0'

   ; Create an "INFO" structure for program info. Store in TLB.

info = { start:start, $         ; ID of START button.
         stop:stop, $           ; ID of STOP button.
         quit:quit, $           ; ID of QUIT button.
         wid:wid, $             ; Graphics window index number of draw widget.
         data:data, $           ; The data to animate.
         delay:0.25, $          ; The amount of time (seconds) to delay between frames.
         stopflag:0, $          ; A stop animation flag.
         nframes:dataSize(3), $ ; The total number of animation frames.
         thisFrame:0 }          ; The current animation frame.

WIDGET_CONTROL, tlb, SET_UVALUE=info, /NO_COPY

   ; Run the program.

XMANAGER, 'xmovie', tlb, GROUP_LEADER=group, Event_Handler='XMOVIE_EVENT'
END
