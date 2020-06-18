;+
; NAME:
;       PROCESS
;
; PURPOSE:
;       The purpose of this routine is to demonstrate a simple
;       image processing program that runs as a "color aware"
;       application. The program works on 8-bit, 16-bit, and
;       24-bit displays. The image is displayed in a resizeable
;       graphics window and the window contents can be saved as
;       GIF or JPEG images. The application can "protect" its
;       colors from other applications that change colors. The
;       color protection is implemented via draw widget EXPOSE
;       events or the Refresh Colors button in the menubar.
;
; AUTHOR:
;       FANNING SOFTWARE CONSULTING
;       David Fanning, Ph.D.
;       2642 Bradbury Court
;       Fort Collins, CO 80521 USA
;       Phone: 970-221-0438
;       E-mail: davidf@dfanning.com
;       Coyote's Guide to IDL Programming: http://www.dfanning.com
;
; CATEGORY:
;       Widgets.
;
; CALLING SEQUENCE:
;       PROCESS
;
; INPUTS:
;       image: An optional 2D image. The image is always scaled.
;
; KEYWORD PARAMETERS:
;       BOTTOM: The lowest color index of the colors to be changed.
;
;       GROUP: The group leader for this program. When the group leader
;       is destroyed, this program will be destroyed.
;
;       NCOLORS: This is the number of colors to load when a color table
;       is selected.
;
; COMMON BLOCKS:
;       None.
;
; SIDE EFFECTS:
;       This is a non-blocking widget.
;
; RESTRICTIONS:
;       The GETIMAGE and XCOLORS programs from the Fanning Consulting
;       web page are required: http://www.dfanning.com/.
;
; EXAMPLE:
;       To run the program, type:
;
;       PROCESS
;
; MODIFICATION HISTORY:
;       Written by:     David Fanning, 13 April 97.
;       Fixed a bug in the TVImage call. 25 June 97. DWF.
;       Extensively modified to incorporate changing ideas
;         about color protection and operation on 8-bit and
;         24-bit displays. 19 Oct 97. DWF.
;       Whoops. Forgot to make *this* draw widget the current
;         graphics window. 15 Nov 97. DWF.
;       IDL 5.1 changed the way color decomposition works. Had
;         to find a fix for this. 25 May 98. DWF.
;-

PRO Process_Quit, event

   ; This event handler destoys the program's TLB.

Widget_Control, event.top, /Destroy
END ;----------------------------------------------------------------------



PRO Process_Cleanup, tlb

   ; Clean up routine for this program. Delete pixmaps
   ; and free pointers. Come here when the program is
   ; destroyed.

Widget_Control, tlb, Get_UValue=info, /No_Copy
IF N_Elements(info) EQ 0 THEN BEGIN
   Heap_GC
   RETURN
ENDIF

WDelete, info.pixid
Ptr_Free, info.image
END ;----------------------------------------------------------------------



PRO Process_Open_File, event

   ; This event handler is responsible for opening new
   ; image files.

Widget_Control, event.top, Get_UValue=info, /No_Copy

Catch, error
IF error NE 0 THEN BEGIN
   IF error EQ -77 THEN BEGIN
      ok = Dialog_Message('Program GETIMAGE required from Coyote Library.')
      Widget_Control, event.top, Set_UValue=info, /No_Copy
      RETURN
   ENDIF
   ok = Dialog_Message(!Err_String)
   Widget_Control, event.top, Set_UValue=info, /No_Copy
   RETURN
ENDIF

   ; Expose events off until after GetImage is destroyed.

Widget_Control, info.drawID, Draw_Expose_Events=0
image = GetImage(Parent=event.top)
Widget_Control, info.drawID, Draw_Expose_Events=1

s = Size(image)
IF s[0] EQ 2 THEN BEGIN
   Ptr_Free, info.image
   info.image = Ptr_New( BytScl(image, Top=info.ncolors-1) + info.bottom )
ENDIF

IF s[0] GT 2 THEN $
   ok = Dialog_Message('Image parameter must be 2D.')

   ; Display in window as last processed.

pseudoEvent = {WIDGET_BUTTON, ID:info.lastid, TOP:event.top, $
      HANDLER:0L, SELECT:1}
Widget_Control, info.lastid, Send_Event=pseudoEvent

Widget_Control, event.top, Set_UValue=info, /No_Copy
END ;----------------------------------------------------------------------



PRO Process_GIF, event

   ; This event handler is responsible for creating a GIF file.

thisFile = Dialog_Pickfile(/Write, File='idl.gif')
IF thisFile EQ '' THEN RETURN

Widget_Control, event.top, Get_UValue=info, /No_Copy
WSet, info.wid
thisImage = TVRD()
TVLCT, r, g, b, /Get
Write_GIF, thisFile, thisImage, r, g, b
Widget_Control, event.top, Set_UValue=info, /No_Copy
END ;----------------------------------------------------------------------



PRO Process_JPEG, event

   ; This event handler is responsible for creating a JPEG file.

thisFile = Dialog_Pickfile(/Write, File='idl.jpg')
IF thisFile EQ '' THEN RETURN

Widget_Control, event.top, Get_UValue=info, /No_Copy
WSet, info.wid
thisImage = TVRD()
TVLCT, r, g, b, /Get

   ; Create 24-bit color image.

image24 = BytArr(3, info.xsize, info.ysize)
image24(0,*,*) = r(thisImage)
image24(1,*,*) = g(thisImage)
image24(2,*,*) = b(thisImage)

   ; Write file. Normal image quality in this lossy format.

Write_JPEG, thisFile, image24, True=1, Quality=75
Widget_Control, event.top, Set_UValue=info, /No_Copy
END ;----------------------------------------------------------------------



PRO Process_Refresh_Colors, event

   ; This event handler is responsible for refreshing
   ; colors in the event they do not refresh from EXPOSE events.

Widget_Control, event.top, Get_UValue=info, /No_Copy

WSet, info.wid
TVLCT, info.r, info.g, info.b, info.bottom
IF info.true THEN BEGIN
   Device, Decomposed=0
   pseudoEvent = {WIDGET_BUTTON, ID:info.lastid, TOP:event.top, $
      HANDLER:0L, SELECT:1}
   Widget_Control, info.lastid, Send_Event=pseudoEvent
ENDIF
Widget_Control, event.top, Set_UValue=info, /No_Copy
END ;----------------------------------------------------------------------



PRO Process_Draw_Events, event

   ; This event handler is responsible for color processing
   ; and EXPOSE events.

Widget_Control, event.top, Get_UValue=info, /No_Copy

WSet, info.wid
Device, Decomposed=0

thisEvent = Tag_Names(event, /Structure)

   ; Reload colors for this widget. This is necessary
   ; if, for example, colors are changed from the command line.
   ; If this is a 24-bit display, the graphic must also be
   ; redisplayed by generating a button event.

IF thisEvent EQ 'XCOLORS_LOAD' THEN BEGIN
   WSET, info.wid
   info.r = event.r(info.bottom:info.ncolors-1+info.bottom)
   info.g = event.g(info.bottom:info.ncolors-1+info.bottom)
   info.b = event.b(info.bottom:info.ncolors-1+info.bottom)
   IF info.true THEN BEGIN
      Device, Decomposed=0
      pseudoEvent = {WIDGET_BUTTON, ID:info.lastid, TOP:event.top, $
         HANDLER:0L, SELECT:1}
      Widget_Control, info.lastid, Send_Event=pseudoEvent
   ENDIF
ENDIF

   ; Expose events processed here.

IF thisEvent EQ 'WIDGET_DRAW' THEN BEGIN ; Button events

   TVLCT, info.r, info.g, info.b, info.bottom

   IF event.type EQ 4 THEN BEGIN ; Expose events only

      IF NOT info.true THEN BEGIN ; 8-bit
         WSet, info.wid
         DEVICE, COPY=[0, 0,info.xsize, info.ysize, 0, 0, info.pixid]

      ENDIF ELSE BEGIN ; 24-bit
         Device, Decomposed=0
         pseudoEvent = {WIDGET_BUTTON, ID:info.lastid, TOP:event.top, $
            HANDLER:0L, SELECT:1}
         Widget_Control, info.lastid, Send_Event=pseudoEvent
      ENDELSE

   ENDIF
ENDIF

Widget_Control, event.top, Set_UValue=info, /No_Copy
END ;----------------------------------------------------------------------



PRO Process_Resize, event

   ; This event handler is responsible for resize events.

Widget_Control, event.top, Get_UValue=info, /No_Copy

   ; Resize the draw widget. Update info sizes.

Widget_Control, info.drawID, XSize=event.x, YSize=event.y
info.xsize = event.x
info.ysize = event.y

   ; Delete and recreate the pixmap if necessary.

IF NOT info.true THEN BEGIN
   WDelete, info.pixid
   Window, /Pixmap, /Free, XSize=event.x, YSize=event.y
   info.pixid = !D.Window
ENDIF

   ; Execute the last command by building an image processing
   ; event to send to the Image Processing event handler.

pseudoEvent = {WIDGET_BUTTON, ID:info.lastid, TOP:event.top, $
   HANDLER:0L, SELECT:1}
Widget_Control, info.lastid, Send_Event=pseudoEvent

Widget_Control, event.top, Set_UValue=info, /No_Copy
END ;----------------------------------------------------------------------



PRO Process_Load_Colors, event

   ; This event handler responds to the Data Colors button
   ; by calling XCOLORS. A "callback" ID is furnished with the
   ; NOTIFYID keyword.

Catch, error
IF error NE 0 THEN BEGIN
   IF error EQ -77 THEN BEGIN
      ok = Dialog_Message('Program XCOLORS required from Coyote Library.')
      Widget_Control, event.top, Set_UValue=info, /No_Copy
      RETURN
   ENDIF
   ok = Dialog_Message(!Err_String)
   Widget_Control, event.top, Set_UValue=info, /No_Copy
   RETURN
ENDIF

Widget_Control, event.top, Get_UValue=info, /No_Copy
XColors, Group=event.top, NColors=info.ncolors, Title='Process Colors', $
   Bottom=info.bottom, NotifyID=[info.drawID, event.top]
Widget_Control, event.top, Set_UValue=info, /No_Copy
END ;----------------------------------------------------------------------



PRO Process_Image_Processing, event

   ; This event handler responds to Image Processing button
   ; events.

Widget_Control, event.top, Get_UValue=info, /No_Copy

   ; Make the pixmap active if in 8-bit. Otherwise,
   ; the display window is the active window.

IF NOT info.true THEN WSet, info.pixid ELSE BEGIN
    WSet, info.wid
    Device, Decomposed=0
ENDELSE

   ; Button events handled here. Branch on button value.

xsize = info.xsize
ysize = info.ysize
Widget_Control, event.id, Get_Value=buttonValue
CASE buttonValue OF
     'Sobel'   : TV, Sobel(Congrid(*info.image, xsize, ysize, /Interp))
     'Roberts' : TV, Roberts(Congrid(*info.image, xsize, ysize, /Interp))
     'Boxcar'  : TV, Smooth(Congrid(*info.image, xsize, ysize, /Interp), 7)
     'Median'  : TV, Median(Congrid(*info.image, xsize, ysize, /Interp), 7)
     'Original Image': TV, Congrid(*info.image, xsize, ysize, /Interp)
ENDCASE

   ; If 8-bit copy from pixmap to display window.

IF NOT info.true THEN BEGIN
   WSet, info.wid
   DEVICE, Copy=[0, 0, xsize, ysize, 0, 0, info.pixid]
ENDIF

   ; Update last event ID. This is how you know what to
   ; put in the display window when you resize it.

info.lastID = event.id

Widget_Control, event.top, Set_UValue=info, /No_Copy
END ;----------------------------------------------------------------------



PRO PROCESS, image, NColors=ncolors, Bottom=bottom, Group=group

   ; This program demonstrates simple image processing capabilities.

On_Error, 1

   ; Must have IDL 5 because of pointers and other functionality.

thisRelease = StrMid(!Version.Release, 0, 1)
IF thisRelease NE '5' THEN BEGIN
   ok = Widget_Message('This program requires IDL 5 functionality. Sorry.')
   RETURN
ENDIF

   ; Load an image if necessary.

IF N_Params(image) EQ 0 THEN BEGIN
  filename = Filepath(Subdir=['examples','data'], 'worldelv.dat')
  OPENR, lun, filename, /Get_LUN
  image = BytArr(360,360)
  ReadU, lun, image
  FREE_LUN, lun
ENDIF

   ; Is the image 2D?

s = SIZE(image)
IF s[0] NE 2 THEN Message, 'Image parameter must be 2D'
xsize = s[1]
ysize = s[2]

   ; 8-bit color or something else?

Window, XSize=10, YSize=10, /Pixmap, /Free
WDelete, !D.Window
IF !D.N_Colors GE 256 THEN BEGIN
    true = 1
    Device, Decomposed=0
ENDIF ELSE true = 0

   ; Check keywords.

IF N_Elements(ncolors) EQ 0 THEN ncolors = (256 < !D.N_Colors)
IF N_Elements(bottom) EQ 0 THEN bottom = 0

   ; Scale the image and load it into a pointer.

thisImage = Ptr_New(BytScl(image, Top=ncolors-1) + bottom)

   ; Create a top-level base for this program. Give it a title.
   ; Create a menubase to hold menu buttons.

tlb = Widget_Base(Title='Resizeable Image Processing Program', Column=1, $
   /TLB_Size_Events, MBar=menubase)

   ; Create a "File" pull-down menu to hold file operation processes.

filer = Widget_Button(menubase, Menu=1, Value='File')

openit = Widget_Button(filer, Value='Open...', Event_Pro='Process_Open_File')
gif = Widget_Button(filer, Value='Save as GIF file...', $
   Event_Pro='Process_GIF')
gif = Widget_Button(filer, Value='Save as JPEG file...', $
   Event_Pro='Process_JPEG')
quitter = Widget_Button(filer, Value='Quit', /Separator, $
   Event_Pro='Process_Quit')

   ; Create a "Colors" button. The Refresh Button is necessary because
   ; it is not possible to "protect" program colors in every instance.
   ; In cases in which program colors are not protected, this button
   ; will do the job.

colors = Widget_Button(menubase, Menu=2, Value='Colors')
dataColors = Widget_Button(colors, Value='Data Colors', $
   Event_Pro='Process_Load_Colors')
refreshColors = Widget_Button(colors, Value='Refresh Colors', $
   Event_Pro='Process_Refresh_Colors')

   ; Create an "Image Processing" pull-down menu. Separate event
   ; handler for all image processing buttons.

processing = Widget_Button(menubase, Menu=1, Value='Image Processing', $
   Event_Pro='Process_Image_Processing')

   edge = Widget_Button(processing, Menu=1, Value='Edge Enhancement')
      sobel = Widget_Button(edge, Value='Sobel')
      robert = Widget_Button(edge, Value='Roberts')

   smooth = Widget_Button(processing, Menu=1, Value='Image Smoothing')
      boxcar = Widget_Button(smooth, Value='Boxcar')
      median = Widget_Button(smooth, Value='Median')

   orig = Widget_Button(processing, Value='Original Image')

   ; Create a drawbase as a child of the top-level base.

drawbase = Widget_Base(tlb, Row=1)

   ; Create a draw widget to hold the image. Separate event handler.
   ; Retain equals 0 to facility expose events. Expose event is not
   ; set until AFTER graphic is displayed in the draw widget to avoid
   ; having two events upon program startup.

drawID = Widget_Draw(drawbase, XSize=xsize, YSize=ysize, $
   Retain=0, Event_Pro='Process_Draw_Events')

   ; Realize the widget program

Widget_Control, tlb, /Realize

   ; Get the value of the draw widget, which is its window index number.
   ; Can only be done AFTER the draw widget is realized. Also, turn
   ; expose events on.

Widget_Control, drawID, Get_Value=wid, Draw_Expose_Events=1

   ; Make the draw widget the current graphics window

WSet, wid
Device, Decomposed=0

   ; Load a colortable and display the image

LoadCT, 4, NColors=ncolors, Bottom=bottom, /Silent

   ; Get color vectors.

TVLCT, r, g, b, /Get

   ; Display scaled image by sending an event to the
   ; button event handler.

pseudoEvent = {WIDGET_BUTTON, ID:orig, TOP:tlb, $
   HANDLER:0L, SELECT:1}
Widget_Control, orig, Send_Event=pseudoEvent

   ; If 8-bit, create a pixmap for expose event recovery.

IF NOT true THEN BEGIN
   Window, /Pixmap, /Free, XSize=xsize, YSize=ysize
   pixid = !D.Window
   Device, Copy=[0, 0, xsize, ysize, 0, 0, wid]
ENDIF ELSE pixid = -1

   ; Create info structure with program information.

info = {  image:thisImage , $             ; The scaled image data.
          wid:wid, $                      ; The window index number.
          pixid:pixid, $                  ; The pixmap window index number.
          true:true, $                    ; True-color flag.
          r:r(bottom:ncolors-1+bottom), $ ; The red colors.
          g:g(bottom:ncolors-1+bottom), $ ; The green colors
          b:b(bottom:ncolors-1+bottom), $ ; The blue colors
          xsize:xsize, $                  ; Current xsize of draw widget.
          ysize:ysize, $                  ; Current ysize of draw widget.
          ncolors:ncolors, $              ; Length of color vectors.
          bottom:bottom, $                ; Starting color index.
          drawID:drawID, $                ; Draw widget ID.
          lastID:orig}                    ; ID of last command button.

Widget_Control, tlb, Set_UValue=info, /No_Copy

   ; Call XManager to set up the Event Loop and manage the program.

XManager, 'process', tlb, Event_Handler='Process_Resize', $
      /No_Block, Group=group
END ;----------------------------------------------------------------------