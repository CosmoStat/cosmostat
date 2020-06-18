;+
; NAME:
;   WINDOWIMAGE
;
; PURPOSE:
;
;       The purpose of this routine is to demonstrate how to interactively
;       adjust the contrast and brightness (also called the window and level)
;       of an image.
;
; AUTHOR:
;
;   FANNING SOFTWARE CONSULTING
;   David Fanning, Ph.D.
;   1645 Sheely Drive
;   Fort Collins, CO 80526 USA
;   Phone: 970-221-0438
;   E-mail: davidf@dfanning.com
;   Coyote's Guide to IDL Programming: http://www.dfanning.com/
;
; CATEGORY:
;
;       Graphics, Widgets.
;
; CALLING SEQUENCE:
;
;       WINDOWIMAGE, image
;
; INPUTS:
;
;       image:    A 2D image array. Typically, the image will be a 16-bit medical
;                 image of some type. If this parameter is not present, a 16-bit
;                 dicom image from the examples/data directory will be selected.
;
; KEYWORD PARAMETERS:
;
;       COLORTABLE:    The index of an IDL-supplied color table to load. Gray-scale by default.
;
;
; COMMON BLOCKS:
;
;       None.
;
; SIDE EFFECTS:
;
;       None.
;
; RESTRICTIONS:
;
;       The TVIMAGE and COLORBAR__DEFINE programs from the Coyote Library are
;       required to run this program:
;
;          http://www.dfanning.com/programs/tvimage.pro
;          http://www.dfanning.com/programs/colorbar__define.pro
;
; EXAMPLE:
;
;       To display a 16-bit medical image, type:
;
;       IDL> WINDOWIMAGE
;
; MODIFICATION HISTORY:
;
;       Written by: David Fanning, 7 November 2001.
;       Modified the windowing algorithm to keep the window within
;          the data range. Added colorbar and pixmap for smoother
;          updating. 18 Nov 2001. DWF.
;       Fixed a small problem with object cleanup. Added window and
;          level readouts, window resizing. 23 Nov 2001. DWF.
;-
;
;###########################################################################
;
; LICENSE
;
; This software is OSI Certified Open Source Software.
; OSI Certified is a certification mark of the Open Source Initiative.
;
; Copyright © 2001 Fanning Software Consulting.
;
; This software is provided "as-is", without any express or
; implied warranty. In no event will the authors be held liable
; for any damages arising from the use of this software.
;
; Permission is granted to anyone to use this software for any
; purpose, including commercial applications, and to alter it and
; redistribute it freely, subject to the following restrictions:
;
; 1. The origin of this software must not be misrepresented; you must
;    not claim you wrote the original software. If you use this software
;    in a product, an acknowledgment in the product documentation
;    would be appreciated, but is not required.
;
; 2. Altered source versions must be plainly marked as such, and must
;    not be misrepresented as being the original software.
;
; 3. This notice may not be removed or altered from any source distribution.
;
; For more information on Open Source Software, visit the Open Source
; web site: http://www.opensource.org.
;
;###########################################################################


PRO WindowImage_Resize, event

; This event handler resizes the widget.

   ; Get the info structure.

Widget_Control, event.top, Get_UValue=info, /No_Copy

   ; Calculate the new draw widget size.

newx = info.labelxsize > event.x
newy = event.y - info.labelysize

   ; Delete the old pixmap and create a new one

WDelete, info. pixmap
Window, XSize=newx, YSize=newy, /Pixmap, /Free
info.pixmap = !D.Window

   ; Resize the draw widget.

Widget_Control, info.drawID, Draw_XSize=newx, Draw_YSize=newy

   ; Restore original contrast and brightness values.

contrast = info.contrast
brightness = info.brightness
level = (1-brightness/100.)*(info.maxVal - info.minVal) + info.minVal
width = (1-contrast/100.)*(info.maxVal - info.minVal)

   ; Calculate display minimum and maximum.

displayMax = (level + (width / 2))
displayMin = (level - (width / 2))

   ; Adjust the window to keep it within the data range.

IF displayMax GT info.maxval THEN BEGIN
   difference = Abs(displayMax - info.maxval)
   displayMax = displayMax - difference
   displayMin = displayMin - difference
ENDIF
IF displayMin LT info.minval THEN BEGIN
   difference = Abs(info.minval - displayMin)
   displayMin = displayMin + difference
   displayMax = displayMax + difference
ENDIF

   ; Display the image.

WSet, info.pixmap
TVImage, BytScl(*info.image, Min=displayMin, Max=displayMax, TOP=!D.Table_Size-1), $
        Position=[0.1, 0.0, 1.0, 0.8], /Keep_Aspect
info.cbar->Clamp, [displayMin,displayMax], /Draw

   ; Update the label.

windowwidth = Fix([displaymin, displaymax])
windowlevel = Fix(level)
txt = 'Window: [' + StrTrim(windowwidth[0],2) + ',' + StrTrim(windowwidth[1],2) + ']  Level: ' + StrTrim(windowlevel,2)
XYOuts, 0.5, 0.95, /Normal, Alignment=0.5, Font=0, txt

WSet, info.wid
Device, Copy=[0, 0, !D.X_Size, !D.Y_Size, 0, 0, info.pixmap]

Widget_Control, event.top, Set_UValue=info, /No_Copy

END ;---------------------------------------------------------------------------



PRO WindowImage_OriginalSettings, event

; This event handler resotres the original window settings.

   ; Get the info structure.

Widget_Control, event.top, Get_UValue=info, /No_Copy

   ; Restore original contrast and brightness values.

contrast = 25
brightness = 75
level = (1-brightness/100.)*(info.maxVal - info.minVal) + info.minVal
width = (1-contrast/100.)*(info.maxVal - info.minVal)

   ; Calculate display minimum and maximum.

displayMax = (level + (width / 2))
displayMin = (level - (width / 2))

   ; Adjust the window to keep it within the data range.

IF displayMax GT info.maxval THEN BEGIN
   difference = Abs(displayMax - info.maxval)
   displayMax = displayMax - difference
   displayMin = displayMin - difference
ENDIF
IF displayMin LT info.minval THEN BEGIN
   difference = Abs(info.minval - displayMin)
   displayMin = displayMin + difference
   displayMax = displayMax + difference
ENDIF

   ; Display the image.

WSet, info.pixmap
TVImage, BytScl(*info.image, Min=displayMin, Max=displayMax, TOP=!D.Table_Size-1), $
        Position=[0.1, 0.0, 1.0, 0.8], /Keep_Aspect, /NoInterpolation
info.cbar->Clamp, [displayMin,displayMax], /Draw

   ; Update the label.

windowwidth = Fix([displaymin, displaymax])
windowlevel = Fix(level)
txt = 'Window: [' + StrTrim(windowwidth[0],2) + ',' + StrTrim(windowwidth[1],2) + ']  Level: ' + StrTrim(windowlevel,2)
XYOuts, 0.5, 0.95, /Normal, Alignment=0.5, Font=0, txt

WSet, info.wid
Device, Copy=[0, 0, !D.X_Size, !D.Y_Size, 0, 0, info.pixmap]

   ; Update the current contrast and brightness values.

info.contrast = contrast
info.brightness = brightness

Widget_Control, event.top, Set_UValue=info, /No_Copy

END ;---------------------------------------------------------------------------



PRO WindowImage_Quit, event

; Destroy the program.

Widget_Control, event.top, /Destroy
END ;---------------------------------------------------------------------------



PRO WindowImage_Cleanup, tlb

; This is where program clean-up occurs.

Widget_Control, tlb, Get_UValue=info, /No_Copy
IF N_Elements(info) NE 0 THEN BEGIN
   Ptr_Free, info.image
   Obj_Destroy, info.cbar
   WDelete, info.pixmap
ENDIF
END ;---------------------------------------------------------------------------



PRO WindowImage_DrawEvents, event

; This event handler responds to cursor events in the
; graphics window.

   ; Error handing.

Catch, theError
IF theError NE 0 THEN BEGIN
   Catch, /Cancel
   ok = Dialog_Message(!Error_State.Msg)
   IF N_Elements(info) NE 0 THEN $
      Widget_Control, event.top, Set_UValue=info, /No_Copy
   RETURN
ENDIF

   ; Get the info structure.

Widget_Control, event.top, Get_UValue=info, /No_Copy

   ; What kind of event is this?

possibleEvent = ['DOWN', 'UP', 'MOTION', 'VIEWPORT', 'EXPOSE']
thisEvent = possibleEvent[event.type]

   ; Branch on type of event.

CASE thisEvent OF

   'DOWN': BEGIN

         ; Set the initial (x,y) point. Turn motion events on.

      info.x = event.x
      info.y = event.y
      Widget_Control, info.drawID, Draw_Motion_Events=1
      END

   'UP': BEGIN

         ; Turn motion events off. Clear any motion events that might have queued.

      Widget_Control, info.drawID, Draw_Motion_Events=0
      Widget_Control, info.drawID, /Clear_Events

         ; Calculate new contrast, brightness, level, and width parameters.

      contrast = 0 > ((info.y - event.y) * info.cstep + info.contrast) < 99
      brightness = 0 > ((info.x - event.x) * info.bstep + info.brightness) < 100
      level = (1-brightness/100.)*(info.maxVal - info.minVal) + info.minVal
      width = (1-contrast/100.)*(info.maxVal - info.minVal)

         ; Calculate display minimum and maximum.

      displayMax = (level + (width / 2))
      displayMin = (level - (width / 2))

         ; Adjust the window to keep it within the data range.

      IF displayMax GT info.maxval THEN BEGIN
         difference = Abs(displayMax - info.maxval)
         displayMax = displayMax - difference
         displayMin = displayMin - difference
      ENDIF
      IF displayMin LT info.minval THEN BEGIN
         difference = Abs(info.minval - displayMin)
         displayMin = displayMin + difference
         displayMax = displayMax + difference
      ENDIF

         ; Display the image.

      WSet, info.pixmap
      TVImage, BytScl(*info.image, Min=displayMin, Max=displayMax, TOP=!D.Table_Size-1), $
              Position=[0.1, 0.0, 1.0, 0.8], /Keep_Aspect, /NoInterpolation, /Erase
      info.cbar->Clamp, [displayMin,displayMax], /Draw

         ; Update the label.

      windowwidth = Fix([displaymin, displaymax])
      windowlevel = Fix(level)
      txt = 'Window: [' + StrTrim(windowwidth[0],2) + ',' + StrTrim(windowwidth[1],2) + ']  Level: ' + StrTrim(windowlevel,2)
      XYOuts, 0.5, 0.95, /Normal, Alignment=0.5, Font=0, txt

      WSet, info.wid
      Device, Copy=[0, 0, !D.X_Size, !D.Y_Size, 0, 0, info.pixmap]

         ; Update the current contrast and brightness values.

      info.contrast = contrast
      info.brightness = brightness
      END

   'MOTION': BEGIN

         ; Calculate new contrast, brightness, level, and width parameters.

      contrast = 0 > ((info.y - event.y) * info.cstep + info.contrast) < 99
      brightness = 0 > ((info.x - event.x) * info.bstep + info.brightness) < 100
      level = (1-brightness/100.)*(info.maxVal - info.minVal) + info.minVal
      width = (1-contrast/100.)*(info.maxVal - info.minVal)

         ; Calculate display minimum and maximum.

      displayMax = (level + (width / 2))
      displayMin = (level - (width / 2))

         ; Adjust the window to keep it within the data range.

      IF displayMax GT info.maxval THEN BEGIN
         difference = Abs(displayMax - info.maxval)
         displayMax = displayMax - difference
         displayMin = displayMin - difference
      ENDIF
      IF displayMin LT info.minval THEN BEGIN
         difference = Abs(info.minval - displayMin)
         displayMin = displayMin + difference
         displayMax = displayMax + difference
      ENDIF

         ; Display the image.

      WSet, info.pixmap
      TVImage, BytScl(*info.image, Min=displayMin, Max=displayMax, TOP=!D.Table_Size-1), $
              Position=[0.1, 0.0, 1.0, 0.8], /Keep_Aspect, /NoInterpolation, /Erase
      info.cbar->Clamp, [displayMin,displayMax], /Draw

         ; Update the label.

      windowwidth = Fix([displaymin, displaymax])
      windowlevel = Fix(level)
      txt = 'Window: [' + StrTrim(windowwidth[0],2) + ',' + StrTrim(windowwidth[1],2) + ']  Level: ' + StrTrim(windowlevel,2)
      XYOuts, 0.5, 0.95, /Normal, Alignment=0.5, Font=0, txt

      WSet, info.wid
      Device, Copy=[0, 0, !D.X_Size, !D.Y_Size, 0, 0, info.pixmap]
      END

   ELSE:

ENDCASE

   ; Return the info structure to storage.

Widget_Control, event.top, Set_UValue=info, /No_Copy
END ;---------------------------------------------------------------------------



PRO WindowImage, image, Colortable=colortable

; The purpose of this procedure is demonstrate one way to
; window or level an image. The image argument must be a
; 2D image.

   ; Error handing.

Catch, theError
IF theError NE 0 THEN BEGIN
   Catch, /Cancel
   ok = Dialog_Message(!Error_State.Msg)
   RETURN
ENDIF

   ; Supply an image if one is not provided.

IF N_Elements(image) EQ 0 THEN BEGIN
   filename = Filepath(Subdir=['examples','data'], 'mr_knee.dcm')
   image = Read_Dicom(filename)
   image = Rebin(Temporary(image), 512, 512)
ENDIF

   ; Load a color table.

IF N_Elements(colortable) EQ 0 THEN colortable = 0 ELSE colortable = 0 > colortable < 40
LoadCT, colortable, /Silent

   ; Make sure image is 2D.

dim = Size(Image, /N_Dimensions)
IF dim NE 2 THEN Message, 'Image argument must be 2D image.'

   ; Find size of image.

s = Size(image, /Dimensions)
xsize = s[0]
ysize = s[1]

   ; Create widgets.

tlb = Widget_Base(Title='Image Window/Leveling Example', Column=1, MBar=menuID, Base_Align_Center=1, $
   TLB_Size_Events=1)
fileID = Widget_Button(menuID, Value='File')
button = Widget_Button(fileID, Value='Restore Original Settings', Event_Pro='WindowImage_OriginalSettings')
quitID = Widget_Button(fileID, Value='Quit', Event_Pro="WindowImage_Quit", /Separator)
drawID = Widget_Draw(tlb, XSize=xsize*1.2, YSize=ysize, Button_Events=1, $
   Event_Pro='WindowImage_DrawEvents')
labelID = Widget_Label(tlb, Value='Drag cursor horizontally for BRIGHTNESS, vertically for CONTRAST.', $
   /Dynamic_Resize)

   ; Find the window index number of the graphics window.

Widget_Control, tlb, /Realize
Widget_Control, drawID, Get_Value=wid
WSet, wid

   ; Set up initial parameters. Contrast and brightness values
   ; go from 0 to 100. Start with 25% contrast and 75% brightness.

maxval = Max(image, Min=minVal)
contrast = 25
brightness = 75

   ; Calculate window level and width from contrast/brightness values.

level = (1-brightness/100.)*(maxVal - minVal) + minVal
width = (1-contrast/100.)*(maxVal - minVal)
displayMax = level + (width / 2)
displayMin = level - (width / 2)

   ; Adjust the window to keep it within the data range.

IF displayMax GT maxval THEN BEGIN
   difference = Abs(displayMax - maxval)
   displayMax = displayMax - difference
   displayMin = displayMin - difference
ENDIF
IF displayMin LT minval THEN BEGIN
   difference = Abs(minval - displayMin)
   displayMin = displayMin + difference
   displayMax = displayMax + difference
ENDIF

   ; Create a pixmap window for smoother updating as the cursor
   ; is moved in the window.

Window, /Free, XSize=xsize*1.2, YSize=ysize, /Pixmap
pixmap = !D.Window
TVImage, BytScl(image, Min=displayMin, Max=displayMax, TOP=!D.Table_Size-1), $
   Position=[0.1, 0.0, 1.0, 0.8], /Keep_Aspect, /NoInterpolation
cbar = Obj_New('Colorbar', Range=[minval, maxval], Font=0, Position=[0.15, 0.88, 0.85, 0.93])
cbar->Loadct, colortable
cbar->Draw

         ; Add a label.

windowwidth = Fix([displaymin, displaymax])
windowlevel = Fix(level)
txt = 'Window: [' + StrTrim(windowwidth[0],2) + ',' + StrTrim(windowwidth[1],2) + ']  Level: ' + StrTrim(windowlevel,2)
XYOuts, 0.5, 0.95, /Normal, Alignment=0.5, Font=0, txt

   ; Display the initial image in the display window

WSet, wid
Device, Copy=[0, 0, !D.X_Size, !D.Y_Size, 0, 0, pixmap]

   ; Find the Y screen sizes of the label widget.

g = Widget_Info(labelID, /Geometry)
labelysize = g.scr_ysize
labelxsize = g.scr_xsize

   ; Create an info structure to hold program information.

info = { wid:wid, $                ; Window index number.
         contrast:contrast, $      ; Current contrast value.
         brightness:brightness, $  ; Current brightness value.
         drawid:drawid, $          ; Draw widget identifier.
         pixmap:pixmap, $          ; The pixmap window
         labelysize:labelysize, $  ; The label Y size.
         labelxsize:labelxsize, $  ; The label X size (minimum size for resizing).
         x:-1, $                   ; Initial X value when windowing starts.
         y:-1, $                   ; Initial Y value when windowing starts.
         xsize:xsize, $            ; X size of image.
         ysize:ysize, $            ; Y size of image.
         bstep:xsize/1024., $      ; The amount of brighness change for one pixel movement.
         cstep:ysize/1024., $      ; The amount of contrast change for one pixel movement.
         image:Ptr_New(image), $   ; A pointer to the image data.
         cbar:cbar, $              ; The colorbar object.
         maxval:maxval, $          ; The maximum value of the image.
         minval:minval }           ; The minimum value of the image.

   ; Store the info structure.

Widget_Control, tlb, Set_UValue=info, /No_Copy

XManager, 'windowimage', tlb, /No_Block, Cleanup='WindowImage_Cleanup', $
   Event_Handler='WindowImage_Resize'
END ;---------------------------------------------------------------------------