;+
; NAME:
;       ZOOMBOX
;
; PURPOSE:
;
;       The purpose of this program is to display an image which
;       can be zoomed by drawing a rubberband box on top of it. The
;      "zoomed" image appears in its own window. The program is written
;       in object graphics. Either 8-bit or 24-bit images may be used.
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
;      Image Processing, Widgets.
;
; CALLING SEQUENCE:
;
;      ZOOMBOX, image
;
; INPUTS:
;
;       image:    An 8-bit or 24-bit image.
;
;       GROUP_LEADER: This keyword is used to assign a group leader to this
;                 program. This program will be destroyed when the group
;                 leader is destroyed. Use this keyword if you are calling
;                 ZOOMBOX from another widget program.
;
;       HARDWARE_RENDERING: Set this keyword if you want to render the scene using
;                 the hardware renderer. This is NOT recommended, since this makes
;                 rendering exceedingly slow on most machines.
;
;       INTERPOLATE: Set this keyword to use bilinear interpolation on the sub-image
;                 defined by the box when zooming. The default is to use nearest neighbor
;                 interpolation.
;
; OUTPUTS:
;       None.
;
; COMMON BLOCKS:
;       None.
;
; SIDE EFFECTS:
;
; DEPENDENCIES:
;
;       Uses XCOLORS and PICKCOLOR from the Coyote Library:
;
;          http://www.dfanning.com/programs/xcolors.pro
;          http://www.dfanning.com/programs/pickcolor.pro
;
; PROCEDURE:
;
;       Clicking the left mouse button allows you to drag a rubberband box
;       over the portion of the window you want to zoom into.
;
; EXAMPLE:
;
;        IDL> ZOOMBOX
;
; RESTRICTIONS:
;
;       Object graphics programs are really meant to be used on 24-bit displays.
;       Colors are decidedly ugly in 256 color environments.
;
; MODIFICATION HISTORY:
;
;       Written by: David W. Fanning, 24 April 2000.
;       Modified draw widget event handling.  22 October 2001. DWF.
;-
;###########################################################################
;
; LICENSE
;
; This software is OSI Certified Open Source Software.
; OSI Certified is a certification mark of the Open Source Initiative.
;
; Copyright © 1999-2001 Fanning Software Consulting
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


FUNCTION Normalize, range, Position=position

On_Error, 1
IF N_Params() EQ 0 THEN Message, 'Please pass range vector as argument.'

IF (N_Elements(position) EQ 0) THEN position = [0.0, 1.0] ELSE $
    position=Float(position)
range = Float(range)

scale = [((position[0]*range[1])-(position[1]*range[0])) / $
    (range[1]-range[0]), (position[1]-position[0])/(range[1]-range[0])]

RETURN, scale
END ;-------------------------------------------------------------------------



PRO ZOOMBOX_BOX_COLOR, event
Widget_Control, event.top, Get_UValue=info, /No_Copy

   ; Get the new box color.

newColor = PickColor(Group_Leader=event.top, Cancel=cancelled)

   ; Store it.

IF NOT cancelled THEN info.boxColor = Reform(newColor, 3)
Widget_Control, event.top, Set_UValue=info, /No_Copy
END ;-------------------------------------------------------------------------



PRO ZOOMBOX_COLORS, event

Widget_Control, event.top, Get_UValue=info, /No_Copy

    ; What kind of event is this?

thisEvent = Tag_Names(event, /Structure)
CASE thisEvent OF

    'WIDGET_BUTTON': BEGIN

            ; Load colors so XCOLORS starts up with the correct colors.

        TVLCT, info.r, info.g, info.b
        XColors, Group=event.top, NotifyID=[event.id, event.top], $
            Title='ZoomBox Colors (' + StrTrim(info.drawID,2) + ')'

        ENDCASE

    'XCOLORS_LOAD':BEGIN

            ; Store the new color table vectors from XCOLORS.

        info.r = event.r
        info.g = event.g
        info.b = event.b
        IF Obj_Valid(info.thePalette) THEN $
            info.thePalette->SetProperty, Red=event.r, Green=event.g, Blue=event.b
        info.theWindow->Draw, info.theView

            ; Is a zoom window open? If so, redisplay it as well.

        IF Obj_Valid(info.zoomWindow) THEN info.zoomWindow->Draw, info.zoomView

        ENDCASE

ENDCASE
Widget_Control, event.top, Set_UValue=info, /No_Copy
END ;-------------------------------------------------------------------------



PRO ZOOMBOX_QUITTER, event
Widget_Control, event.top, /Destroy
END ;-------------------------------------------------------------------------



PRO ZOOMBOX_CLEANUP, tlb

   ; The purpose of this program is to clean up the objects created for the
   ; main image window.

Widget_Control, tlb, Get_UValue=info, /No_Copy
IF N_Elements(info) NE 0 THEN BEGIN
    Obj_Destroy, info.theContainer
ENDIF
END ;-------------------------------------------------------------------------



PRO ZOOMBOX_ZOOM_CLEANUP, tlb

   ; The purpose of this program is to is to clean up the objects created for the
   ; zoom window.

Widget_Control, tlb, Get_UValue=zoomContainer
Obj_Destroy, zoomContainer
END ;-------------------------------------------------------------------------



PRO ZOOMBOX_FACTOR, event

   ; The purpose of this event handler is to set the zoom factor.

Widget_Control, event.top, Get_UValue=info, /No_Copy
Widget_Control, event.id, Get_UValue=factor
info.zoomfactor = factor
Widget_Control, event.top, Set_UValue=info, /No_Copy
END ;-------------------------------------------------------------------------



PRO ZOOMBOX_DRAW_EVENTS, event

   ; This event handler responds to events generated by the draw widget.

   ; Get the info structure out of the top-level base.

Widget_Control, event.top, Get_UValue=info, /No_Copy

   ; What type of an event is this?

possibleEventTypes = [ 'DOWN', 'UP', 'MOTION', 'SCROLL', 'EXPOSE' ]
thisEvent = possibleEventTypes(event.type)

CASE thisEvent OF

   'EXPOSE': BEGIN ; Redraw the view.

      info.theWindow->Draw, info.theView

      ENDCASE

   'DOWN': BEGIN

         ; Set the static corners of the box to current
         ; cursor location.

      info.xs = event.x
      info.ys = event.y

         ; Change the event handler for the draw widget and turn MOTION
         ; events ON.

      Widget_Control, event.id, Draw_Motion_Events=1

         ; Initialize the polygon object.

      info.theBox = Obj_New('IDLgrPolyline', Replicate(info.xs, 5) / info.xsize, $
         Replicate(info.ys, 5) / info.ysize, Color=info.boxColor)
      info.boxModel = Obj_New('IDLgrModel')
      info.boxModel->Add, info.theBox

         ; Create an instance of the current view.

      info.theWindow->Draw, info.theView, /Create_Instance
      info.theModel->SetProperty, Hide=0
      info.theView->SetProperty, Transparent=0

         ; Add the boxModel to the view so it can be displayed.

      info.theView->Add, info.boxModel

      ENDCASE

   'UP': BEGIN

         ; If this is an UP event, you need to erase the zoombox, turn motion events
         ; OFF, and draw the "zoomed" plot.

         ; Turn motion events off..

      Widget_Control, event.id, Draw_Motion_Events=0

         ; Draw the "zoomed" image. Start by getting the LAST zoom
         ; box outline. These are indices into image array.

      event.x = 0 > event.x < (info.xsize-1)
      event.y = 0 > event.y < (info.ysize-1)
      x = [info.xs, event.x]
      y = [info.ys, event.y]

         ; Make sure the user didn't just click in the window.

      IF info.xs EQ event.x OR info.ys EQ event.y THEN BEGIN
         Widget_Control, event.top, Set_UValue=info, /No_Copy
         RETURN
      ENDIF

         ; Make sure the x and y values are ordered as [min, max].

      IF info.xs GT event.x THEN x = [event.x, info.xs]
      IF info.ys GT event.y THEN y = [event.y, info.ys]

         ; Set the zoom factor and determine the new X and Y
         ; sizes of the Zoom Window.

      zoomXSize = (x[1] - x[0] + 1) * info.zoomFactor
      zoomYSize = (y[1] - y[0] + 1) * info.zoomFactor

         ; Subset the image, and apply the zoom factor to it.

      info.theImage->GetProperty, Data=image

      CASE info.true OF

         0: BEGIN
            imageSubset = image(x[0]:x[1], y[0]:y[1])
            zoomedImage = Congrid(imageSubset, zoomXSize, zoomYSize, Interp=info.interp)
            ENDCASE

         1: BEGIN
            zoomedImage = BytArr(3, zoomXsize, zoomYSize)
            FOR j = 0,2 DO BEGIN
               slice = image[j,*,*]
               slice = Reform(slice, info.xsize, info.ysize)
               subimage = slice[x[0]:x[1], y[0]:y[1]]
               zoomedImage[j,*,*] = Congrid(subimage, zoomXSize, zoomYSize, Interp=info.interp)
            ENDFOR
            ENDCASE
         2: BEGIN
            zoomedImage = BytArr(zoomXsize, 3, zoomYSize)
            FOR j = 0,2 DO BEGIN
               slice = image[*,j,*]
               slice = Reform(slice, info.xsize, info.ysize)
               subimage = slice[x[0]:x[1], y[0]:y[1]]
               zoomedImage[*,j,*] = Congrid(subimage, zoomXSize, zoomYSize, Interp=info.interp)
            ENDFOR
            ENDCASE
         3: BEGIN
            zoomedImage = BytArr(zoomXsize, zoomYSize, 3)
            FOR j = 0,2 DO BEGIN
               slice = image[*,*,j]
               slice = Reform(slice, info.xsize, info.ysize)
               subimage = slice[x[0]:x[1], y[0]:y[1]]
               zoomedImage[*,*,j] = Congrid(subimage, zoomXSize, zoomYSize, Interp=info.interp)
            ENDFOR
            ENDCASE
      ENDCASE

         ; Create or re-populate (if it already exists) the zoomed image object.

      IF Obj_Valid(info.zoomImage) THEN info.zoomImage->SetProperty, Data=zoomedImage ELSE $
         info.zoomImage = Obj_New('IDLgrImage', zoomedImage, Palette=info.thePalette)

         ; Scale the image into the view.

      info.zoomImage->SetProperty, XCoord_Conv=Normalize([0,zoomXSize]), YCoord_Conv=Normalize([0,zoomYSize])

         ; If the Zoom Window exists, make it the proper size and load
         ; the zoomed image into it. If it does not exists, create it.

      IF Obj_Valid(info.zoomWindow) THEN BEGIN

            ; Zoomed window exists. Make it correct size and load image.

        info.zoomWindow->SetProperty, Dimension=[zoomXSize, zoomYSize]

      ENDIF ELSE BEGIN

            ; Get offset positions for the non-existing zoom window.

         Widget_Control, event.top, TLB_Get_Size=sizes, TLB_Get_Offset=offsets
         xpos = sizes[0] + offsets[0] + 20
         ypos = offsets[1] + 40

            ; Create zoomed window.

         zoomtlb = Widget_Base(Title='Zoomed Image', Group=event.top, Kill_Notify='ZOOMBOX_ZOOM_CLEANUP', $
            XOffset=xpos, YOffset=ypos)
         zoomdraw = Widget_Draw(zoomtlb, XSize=zoomXSize, YSize=zoomYSize, $
            Graphics_Level=2, Retain=2, Renderer=1)
         Widget_Control, zoomtlb, /Realize
         Widget_Control, zoomdraw, Get_Value=zoomWindow
         info.zoomWindow = zoomWindow
         info.zoomView = Obj_New('IDLgrView', Viewplane_Rect=[0,0,1,1])
         info.zoomModel = Obj_New('IDLgrModel')
         info.zoomView->Add, info.zoomModel
         info.zoomModel->Add, info.zoomImage

            ; Create container for cleaning up later.

         zoomContainer = Obj_New('IDL_Container')
         zoomContainer->Add, info.zoomWindow
         zoomContainer->Add, info.zoomModel
         zoomContainer->Add, info.zoomView
         zoomContainer->Add, info.zoomImage
         Widget_Control, zoomtlb, Set_UValue=zoomContainer

      ENDELSE

         ; Update the zoomed image data and draw it.

      info.zoomImage->SetProperty, Data=zoomedImage
      info.zoomWindow->Draw, info.zoomView

         ; Destroy the box object and model.

      Obj_Destroy, info.boxModel
      Obj_Destroy, info.theBox

         ; Set the instancing variables back to default values.

      info.theView->SetProperty, Transparent=0
      info.theModel->SetProperty, Hide=0

         ; Clear any motion events that may have occurred.

      Widget_Control, event.id, /Clear_Events

   ENDCASE

   'MOTION': BEGIN

      ; Erase the old zoom box.

    info.theWindow->Draw, info.theView, /Draw_Instance

      ; Get the dynamic corner of the box.

    info.xd = 0 > event.x < (info.xsize-1)
    info.yd = 0 > event.y < (info.ysize-1)

      ; Re-configure the box coordinates.

   theBoxData = FltArr(2,5)
   theBoxData[0,*] = [info.xs, info.xd, info.xd, info.xs, info.xs] / Float(info.xsize)
   theBoxData[1,*] = [info.ys, info.ys, info.yd, info.yd, info.ys] / Float(info.ysize)

      ; Draw the new zoom box.

   info.theBox->SetProperty, Data=theBoxData
   info.theWindow->Draw, info.theView, /Draw_Instance

   ENDCASE

   ELSE:

ENDCASE

   ; Put the info structure back into its storage location.

Widget_Control, event.top, Set_UValue=info, /No_Copy

END ;-------------------------------------------------------------------------



PRO ZOOMBOX, image, Group_Leader=group, Interpolation=interp, Hardware_Render=hardware

   ; This procedure allows the user to "zoom" into an image
   ; by drawing a box around the part of the image to zoom into.
   ; This is a companion program to ZIMAGE, however this program
   ; is written in object graphics.

   ; On an error condition, return to the main level of IDL.

On_Error, 1

   ; Was an image passed into the procedure?
   ; f not, find one in the Examples directory.

IF N_Params() EQ 0 THEN BEGIN
   filename = Filepath(SubDirectory=['examples','data'], 'worldelv.dat')
   OpenR, lun, filename, /Get_Lun
   image = BytArr(360,360)
   ReadU, lun, image
   Free_Lun, lun
ENDIF

   ; Make sure you have an image of the proper type. 8-bit or 24-bit images are allowed.

nDims = Size(image, /N_Dimensions)
imageSize = Size(image, /Dimensions)
IF nDims EQ 2 THEN BEGIN
   xsize = imageSize[0]
   ysize = imageSize[1]
   true = 0
ENDIF ELSE BEGIN
   IF NDims NE 3 THEN Message, 'Image argument does not appear to be an image. Returning...'
   IF imageSize[0] NE 3 AND imageSize[1] NE 3 AND imageSize[2] NE 3 THEN BEGIN
      Message, 'Image does not appear to be a 24-bit image. Returning...'
   ENDIF
   IF imageSize[0] EQ 3 THEN BEGIN
      xsize = imageSize[1]
      ysize = imageSize[2]
      true = 1
   ENDIF
   IF imageSize[1] EQ 3 THEN BEGIN
      xsize = imageSize[0]
      ysize = imageSize[2]
      true = 2
   ENDIF
   IF imageSize[2] EQ 3 THEN BEGIN
      xsize = imageSize[0]
      ysize = imageSize[1]
      true = 3
   ENDIF
ENDELSE

   ; Do you need a color palette? If so, load gray-scale colortable.

IF true EQ 0 THEN BEGIN
   thePalette = Obj_New('IDLgrPalette')
   thePalette->LoadCT, 0
ENDIF ELSE thePalette = Obj_New()

   ; Create the necessary objects.

theImage = Obj_New('IDLgrImage', image, Palette=thePalette)
theModel = Obj_New('IDLgrModel')
theModel->Add, theImage
theView = Obj_New('IDLgrView', Viewplane_Rect=[0,0,1,1])
theView->Add, theModel
theImage->SetProperty, XCoord_Conv=Normalize([0,xsize]), YCoord_Conv=Normalize([0,ysize])

   ; Store objects in container for later cleaning up.

theContainer = Obj_New('IDL_Container')
theContainer->Add, theModel
theContainer->Add, theView
theContainer->Add, thePalette
theContainer->Add, theImage

   ; Create a top-level base for this program. No resizing of this base.

tlb = Widget_Base(TLB_Frame_Attr=1, MBar=menuID)

   ; Menu items.

controlID = Widget_Button(menuID, Value='Controls')
IF true EQ 0 THEN colors = Widget_Button(controlID, Value='Image Colors...', $
   Event_Pro='ZOOMBOX_COLORS')
boxColorID = Widget_Button(controlID, Value='Zoom Box Color...', $
   Event_Pro='ZOOMBOX_BOX_COLOR')
zoomfactorID = Widget_Button(controlID, Value='Zoom Factor', /Menu, $
   Event_Pro='ZOOMBOX_FACTOR')
button = Widget_Button(zoomfactorID, Value='2X', UValue=2.0)
button = Widget_Button(zoomfactorID, Value='3X', UValue=3.0)
button = Widget_Button(zoomfactorID, Value='5X', UValue=5.0)
button = Widget_Button(zoomfactorID, Value='8X', UValue=8.0)
quitter = Widget_Button(controlID, Value='Quit', $
   Event_Pro='ZOOMBOX_QUITTER', /Separator)

   ; Draw widget. Make sure you use software rendering or rubberband box
   ; drawing can be V-E-R-Y slow.

drawbase = Widget_Base(tlb, Map=1)
drawID = Widget_Draw(drawbase, XSize=xsize, YSize=ysize, Retain=0, $
   Button_Events=1, Event_Pro='ZOOMBOX_DRAW_EVENTS', Graphics_Level=2, $
   Expose_Events=1, Renderer=1-Keyword_Set(hardware))

   ; Realize the program.

Widget_Control, tlb, /Realize

   ; Get the window object.

Widget_Control, drawID, Get_Value=theWindow

   ; Set the title of the window.

Widget_Control, tlb, TLB_Set_Title='Full Size Image (' + StrTrim(drawID,2) + ')'

   ; Get color vectors for this application.

TVLCT, r, g, b, /Get

   ; Select an ARROW cursor.

theWindow->SetCurrentCursor, 'Arrow'

   ; Create an info structure to hold information required by the program.

info = { $
   theImage:theImage, $         ; The original image.
   zoomImage:Obj_New(), $       ; The scaled and resized subimage.
   theView:theView, $           ; The view that will be rendered.
   theModel:theModel, $         ; The image model.
   theContainer:theContainer, $ ; The main container object.
   theWindow:theWindow, $       ; The main image window object.
   thePalette:thePalette, $     ; The palette object.
   theBox:Obj_New(), $          ; The rubberband box object.
   boxColor:[0,255,0], $        ; The box color. Green to start.
   boxModel:Obj_New(), $        ; The box model
   xsize:xsize, $               ; The x size of the image window and of the image.
   ysize:ysize, $               ; The y size of the image window and of the image.
   xs:0, $                      ; X static corner of the zoom box.
   ys:0, $                      ; Y static corner of the zoom box.
   xd:0, $                      ; X dynamic corner of the zoom box.
   yd:0, $                      ; Y dynamic corner of the zoom box.
   r:r, $                       ; The red vector of the color table.
   g:g, $                       ; The green vector of the color table.
   b:b, $                       ; The blue vector of the color table.
   true:true, $                 ; A flag that indicates 24-bit image interleaving.
   drawID:drawID, $             ; The draw widget identifier.
   zoomWindow:Obj_New(), $      ; The zoom window object.
   zoomView:Obj_New(), $        ; The zoom window view object.
   zoomModel:Obj_New(), $       ; The zoom model.
   zoomfactor:2, $              ; The initial zoom factor.
   interp:Keyword_Set(interp) $ ; A flag that indicates interpolation of the subimage.
    }

   ; Store the info structure in the user value of the top-level base.

Widget_Control, tlb, Set_UValue=info, /No_Copy

   ; Register this program and set up the event loop.

XManager, 'zoombox', tlb, Cleanup='ZOOMBOX_CLEANUP', Group_Leader=group, /No_Block
END ;-------------------------------------------------------------------------
