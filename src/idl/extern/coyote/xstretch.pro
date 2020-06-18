;+
; NAME:
;       XSTRETCH
;
; PURPOSE:
;
;       The purpose of this program is to interactively apply a simple
;       linear stretch to an image by moving two lines on a histogram
;       plot of the image. The portion of the image data between the
;       two lines is stretched over the available colors in the color table.
;
; AUTHOR:
;
;       FANNING SOFTWARE CONSULTING
;       David Fanning, Ph.D.
;       1645 Sheely Drive
;       Fort Collins, CO 80526 USA
;       Phone: 970-221-0438
;       E-mail: davidf@dfanning.com
;       Coyote's Guide to IDL Programming: http://www.dfanning.com
;
; CATEGORY:
;
;       Graphics, Widgets
;
; CALLING SEQUENCE:
;
;       XSTRETCH, image
;
; INPUT PARAMETERS:
;
;       image:    The image data to be stretched.It must be 2D. (This may now
;                 be a pointer to the image data rather than the image itself.)
;
; KEYWORD PARAMETERS:
;
;       COLORTABLE: The index of a colortable you would like to load.
;                 The current colortable is used if this keyword is undefined.
;
;       _EXTRA:   This keyword collects any keyword appropriate for the
;                 Plot command.
;
;       GROUP_LEADER: Keyword to assign a group leader (so this program can be
;                 called from within another widget program).
;
;       MAX_VALUE: Keyword to assign a maximun value for the Histogram Plot.
;                 Images with lots of pixels of one color (e.g. black) skew
;                 the histogram. This helps make a better looking plot.
;
;       NCOLORS:  Keyword to assign the number of colors used to display
;                 the image. The default is !D.Table_Size-4.
;
; OUTPUTS:
;       None.
;
; COMMON BLOCKS:
;       None.
;
; SIDE EFFECTS:
;       None.
;
; DEPENDENCIES:
;
;       The following programs are required from the Coyote Library:
;
;           tvimage.pro
;           pswindow.pro
;           psconfig.pro
;           fsc_inputfield.pro
;           fsc_droplist.pro
;           fsc_fileselect.pro
;           fsc_plotwindow.pro
;           fsc_psconfig__define.pro
;
; EXAMPLE:
;
;       If you have a 2D image in the variable "image", you can run this
;       program like this:
;
;       XSTRETCH, image
;
; MODIFICATION HISTORY:
;
;       Written by: David Fanning, April 1996.
;       October, 1996 Fixed a problem with not restoring the color
;          table when the program exited. Substituted a call to XCOLORS
;          instead of XLOADCT.
;       October, 1998. Added NO_BLOCK keyword and modified to work with
;          24-bit color devices.
;       April, 1999. Made lines thicker. Offered default image. DWF.
;       April, 1999. Replaced TV command with TVIMAGE. DWF.
;       April, 1999. Made both windows resizeable. DWF.
;       April, 2000. Made several modifications to histogram plot and to
;          the way colors were handled. Added ability to pass pointer to
;          the image as well as image itself. DWF.
;       February 2001. Removed GIF file support for IDL 5.4 and fixed
;          a problem with cleaning up the pixmap. DWF.
;-
;
;###########################################################################
;
; LICENSE
;
; This software is OSI Certified Open Source Software.
; OSI Certified is a certification mark of the Open Source Initiative.
;
; Copyright © 2000-2001 Fanning Software Consulting.
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

PRO XSTRETCH_IMAGEWINDOWKILLED, imageWindowID

; Turn the Save As, Print, and Image Colors buttons off.

Widget_Control, imageWindowID, Get_UValue=buttonIDs
IF Widget_Info(buttonIDs[0], /Valid_ID) THEN BEGIN
   Widget_Control, buttonIDs[0], Sensitive=0
   Widget_Control, buttonIDs[1], Sensitive=0
   Widget_Control, buttonIDs[2], Sensitive=0
ENDIF
END ;--------------------------------------------------------------------



PRO XSTRETCH_HISTOPLOT, image, Binsize=binsize, Reverse_Indices=r, WID=wid, $
   Color=color, Background=background, Max_Value=maxvalue, _Extra=extra

; This is a utility program to draw a histogram plot in a
; display window.

   ; Catch any error in the histogram display routine.

Catch, theError
;theError=0
IF theError NE 0 THEN BEGIN
   Catch, Cancel=1
   ok = Dialog_Message('Error in XStretch_Histoplot. Returning...')
   Print, ''
   Print, 'XStretch_Histoplot: ' + !Error_State.Msg
   RETURN
ENDIF

   ; Calculate binsize.

range = Max(*image) - Min(*image)
binsize = 1.0 > (range / 255.)

histdata = Histogram(*image, Binsize=binsize, Reverse_Indices=r)

   ; Plot the histogram of the display image.

IF N_Elements(wid) NE 0 THEN WSet, wid
bins = Findgen(N_Elements(histdata)) * binsize + Min(*image)
xrange = [Min(bins), Max(bins)]
yrange = [0,maxValue]
Plot, bins, histdata, YTitle='Pixel Density', $
   Background=background, Color=color, /NoData, $
   XRange=xrange, XStyle=1, Max_Value=maxValue,  $
   XTitle='Image Value', Title='Image Histogram', _Extra=extra, $
   XTickformat='(I6)', YTickformat='(I6)', YRange=yrange, YStyle=1
FOR j=0L,N_Elements(bins)-2 DO BEGIN
   PlotS, [bins[j], bins[j], bins[j+1], bins[j+1]], $
          [0, histdata[j] < !Y.CRange[1], histdata[j] < !Y.CRange[1], 0], Color=color
ENDFOR

END ;--------------------------------------------------------------------------------



PRO XSTRETCH_SAVEAS, event

   ; Errors caused by incorrect IDL versions or missing Coyote files.

Catch, theError
IF theError NE 0 THEN BEGIN
   Catch, /Cancel
   ok = Dialog_Message(['Sorry. This code uses keywords', $
                        'or programs not available in this', $
                        'version of IDL. Returning...'])
   Print, !Error_State.Msg
   IF N_Elements(info) NE 0 THEN Widget_Control, event.top, Set_UValue=info, /No_Copy
   RETURN
ENDIF

   ; Save as various file types.

Widget_Control, event.top, Get_UValue=info, /No_Copy
Widget_Control, event.id, Get_UValue=saveAsType

   ; Preliminary information.

WSet, info.windex
xsize = !D.X_Size
ysize = !D.Y_Size
Device, Get_Visual_Depth=thisDepth, Get_Decomposed=decomposedState

   ; What kind of file do you want?

CASE saveAsType OF

   'JPEG': BEGIN

         filename = Dialog_Pickfile(/Write, File='xstretch.jpg')
         IF filename NE "" THEN BEGIN
            IF thisDepth GT 8 THEN BEGIN
               Device, Decomposed=1
               image24 = TVRD(True=1)
               Write_JPEG, filename, image24, True=1
               Device, Decomposed=decomposedState
            ENDIF ELSE BEGIN
               image24 = BytArr(3, xsize, ysize)
               image2d = TVRD()
               TVLCT, r, g, b, /Get
               image24[0,*,*] = r[image2d]
               image24[1,*,*] = g[image2d]
               image24[2,*,*] = b[image2d]
               Write_JPEG, filename, image24, True=1
            ENDELSE
         ENDIF

         ENDCASE

   'TIFF': BEGIN

         filename = Dialog_Pickfile(/Write, File='xstretch.tif')
         IF filename NE "" THEN BEGIN
            IF thisDepth GT 8 THEN BEGIN
               Device, Decomposed=1
               image24 = TVRD(True=1)
               image24 = Reverse(Temporary(image24),3)
               Write_TIFF, filename, image24, 1
               Device, Decomposed=decomposedState
            ENDIF ELSE BEGIN
               image24 = BytArr(3, xsize, ysize)
               image2d = TVRD()
               TVLCT, r, g, b, /Get
               image24[0,*,*] = r[image2d]
               image24[1,*,*] = g[image2d]
               image24[2,*,*] = b[image2d]
               image24 = Reverse(Temporary(image24),3)
               Write_TIFF, filename, image24, 1
            ENDELSE
         ENDIF

         ENDCASE

   'GIF': BEGIN

         filename = Dialog_Pickfile(/Write, File='xstretch.gif')
         IF filename NE "" THEN BEGIN
            IF thisDepth GT 8 THEN BEGIN
               Device, Decomposed=1
               image24 = TVRD(True=1)
               image2d = Color_Quan(image24, 1, r, g, b)
               Write_GIF, filename, image2d, r, g, b
               Device, Decomposed=decomposedState
            ENDIF ELSE BEGIN
               image2d = TVRD()
               TVLCT, r, g, b, /Get
               Write_GIF, filename, image2d, r, g, b
            ENDELSE
         ENDIF

         ENDCASE

   'PS': BEGIN

         configureIt = PSConfig(Group_Leader=event.top, Cancel=cancelled, Color=1, $
            Filename='xstretch.ps')
         IF NOT cancelled THEN BEGIN
               thisDevice = !D.Name
               Set_Plot, 'PS', /Copy
               Device, _Extra=configureIt
                  displayImage = BytScl(*info.image, Top=info.ncolors-1,  $
                  Max=info.maxThresh, Min=info.minThresh)
               TVImage, displayImage, _Extra=info.extra
               Device, /Close_File
               Set_Plot, thisDevice
         ENDIF

         ENDCASE



ENDCASE


   ; Put the info structure back.

Widget_Control, event.top, Set_UValue=info, /No_Copy
END
;-------------------------------------------------------------------



PRO XSTRETCH_PRINT, event

   ; Printing and printer setup handled here.

Widget_Control, event.top, Get_UValue=info, /No_Copy

   ; Configure printer and print if user OKs.

result = Dialog_PrinterSetup()
IF result EQ 1 THEN BEGIN

      ; I want the output on the page to have the same aspect ratio
      ; as I see in the display window.

   WSet, info.windex
   configurePrinter = PSWindow(/Printer)

      ; Print the image.

   thisDevice = !D.Name
   Set_Plot, 'PRINTER'
   Device, _Extra=configurePrinter
   Widget_Control, Hourglass=1
   displayImage = BytScl(*info.image, Top=info.ncolors-1,  $
      Max=info.maxThresh, Min=info.minThresh)
   TVImage, displayImage, _Extra=info.extra
   Widget_Control, Hourglass=0
   Device, /Close_Document
   Set_Plot, thisDevice

ENDIF

   ; Put the info structure back.

Widget_Control, event.top, Set_UValue=info, /No_Copy
END
;-------------------------------------------------------------------



PRO XSTRETCH_PROCESS_EVENTS, event

   ; This event handler ONLY responds to button down events from the
   ; draw widget. If it gets a DOWN event, it does two things: (1) finds
   ; out which threshold line is to be moved, and (2) changes the
   ; event handler for the draw widget to XSTRETCH_MOVELINE.

possibleEventTypes = [ 'DOWN', 'UP', 'MOTION', 'SCROLL' ]
thisEvent = possibleEventTypes(event.type)
IF thisEvent NE 'DOWN' THEN RETURN

      ; Must be DOWN event to get here, so get info structure.

   Widget_Control, event.top, Get_UValue=info, /No_Copy

      ; Convert the device coordinates to data coordinates.
      ; Have to have scaling factors for conversion.

   Wset, info.histo_wid
   !X.S = info.xs
   !Y.S = info.ys
   coords = Convert_Coord(event.x, event.y, 0, /Device, /To_Data)

      ; Is this event close to a line? If not, ignore it.

      ; Click has to be inside the graph in the y direction.

   IF coords(1) LT info.ymin OR coords(1) GT info.ymax THEN BEGIN
      Widget_Control, event.top, Set_UValue=info, /No_Copy
      RETURN
   ENDIF

      ; How close to either line are you?

    closemin = Abs(info.minthresh - coords(0))
    closemax = Abs(info.maxthresh - coords(0))
    IF closemin LE closemax THEN info.lineby = 'MIN' ELSE info.lineby = 'MAX'

       ; If you are not close to a line, goodbye!

    CASE info.lineby OF
       'MIN': BEGIN
              IF closemin GT info.close THEN BEGIN
                  Widget_Control, event.top, Set_UValue=info, /No_Copy
                  RETURN
              ENDIF
              END

       'MAX': BEGIN
              IF closemax GT info.close THEN BEGIN
                  Widget_Control, event.top, Set_UValue=info, /No_Copy
                  RETURN
              ENDIF
              END
    ENDCASE

    ; Change the event handler for the draw widget and turn MOTION
    ; events ON.

 Widget_Control, event.id, Event_Pro='XSTRETCH_MOVELINE', $
    Draw_Motion_Events=1

   ; Put the info structure back into its storage location.

Widget_Control, event.top, Set_UValue=info, /No_Copy
END ; of XSTRETCH_PROCESS_EVENTS *********************************************



PRO XSTRETCH_MOVELINE, event

   ; This event handler continuously draws and erases a threshold line
   ; until it receives an UP event from the draw widget. Then it turns
   ; draw widget motion events OFF and changes the event handler for the
   ; draw widget back to XSTRETCH_PROCESS_EVENTS.

   ; Get the info structure out of the top-level base.

Widget_Control, event.top, Get_UValue=info, /No_Copy

   ; What type of an event is this?

possibleEventTypes = [ 'DOWN', 'UP', 'MOTION', 'SCROLL' ]
thisEvent = possibleEventTypes(event.type)

IF thisEvent EQ 'UP' THEN BEGIN

      ; If this is an UP event, set the draw widget's event handler back
      ; to XSTRETCH_PROCESS_EVENTS, turn MOTION events OFF, and apply the
      ; new threshold parameters to the image.

      ; Erase the last theshold line drawn.

   WSet, info.histo_wid
   !X.S = info.xs
   !Y.S = info.ys
   Device, Copy = [0, 0, info.pix_xsize, info.pix_ysize, 0, 0, info.pixmap]

      ; Turn motion events off and redirect the events to XSTRETCH_PROCESS_EVENTS.

    Widget_Control, event.id, Draw_Motion_Events=0, $
       Event_Pro='XStretch_Process_Events'

      ; Convert the event device coordinates to data coordinates.

   coord = Convert_Coord(event.x, event.y, /Device, /To_Data)

      ; Make sure the coordinate is between the other line and
      ; still inside the plot.

   CASE info.lineby OF
      'MIN': BEGIN
             coord(0) = coord(0) > (info.xmin + 1)
             coord(0) = coord(0) < (info.maxThresh - 1)
             END
      'MAX': BEGIN
             coord(0) = coord(0) > (info.minThresh + 1)
             coord(0) = coord(0) < (info.xmax - 1)
             END
   ENDCASE

      ; Draw both of the threshold lines again.

   CASE info.lineby OF
      'MIN': BEGIN
             PlotS, [coord(0), coord(0)],[info.ymin, info.ymax], $
                Color=info.minColor, Thick=3
             PlotS, [info.maxThresh, info.maxThresh],  $
                [info.ymin, info.ymax], Color=info.maxColor, Thick=3
             info.minthresh = coord(0)
             END
      'MAX': BEGIN
             PlotS, [coord(0), coord(0)],[info.ymin, info.ymax],  $
                 Color=info.maxColor, Thick=3
             PlotS, [info.minThresh, info.minThresh],  $
                [info.ymin, info.ymax], Color=info.minColor, Thick=3
             info.maxthresh = coord(0)
             END
   ENDCASE

   ; Update the image display by appling the threshold parameters.
   ; Be sure the image draw widget is still around. Make it if it isn't.

IF Widget_Info(info.image_draw, /Valid_ID) THEN BEGIN
   WSet, info.windex
   displayImage = BytScl(*info.image, Top=info.ncolors-1,  $
      Max=info.maxThresh, Min=info.minThresh)
   TVImage, displayImage, _Extra=info.extra

ENDIF ELSE BEGIN

   imageSize = Size(*info.image)
   xsize = imageSize(1)
   ysize = imageSize(2)
   aspect = Float(xsize)/ysize
   IF xsize GT 512 OR ysize GT 512 THEN BEGIN
      IF xsize GT ysize THEN BEGIN
         xsize = 512
         ysize = Fix(xsize / aspect)
      ENDIF ELSE BEGIN
         ysize = 512
         xsize = ysize / aspect
      ENDELSE
   ENDIF
   Widget_Control, event.top, TLB_Get_Offset=offsets, TLB_Get_Size=sizes
   xoff = offsets[0] + sizes[0] + 20
   yoff = offsets[1]
   image_tlb = Widget_Base(Row=1, Group=event.top, Title='XStretch Image', $
      XOffSet=xoff, YOffSet=yoff, TLB_Size_Events=1)
   image_draw = Widget_Draw(image_tlb, XSize=xsize, YSize=ysize)
   Widget_Control, image_tlb, /Realize, Set_UValue=event.top
   Widget_Control, image_draw, Get_Value=windex
   info.image_draw = image_draw
   info.windex = windex
   displayImage = BytScl(*info.image, Top=info.ncolors-1,  $
      Max=info.maxThresh, Min=info.minThresh)
   TVImage, displayImage, _Extra=info.extra
   XManager, 'xstretch_image', image_tlb, Event_Handler='XStretch_Image_Resize', /No_Block
   Widget_Control, info.saveas, Sensitive=1
   Widget_Control, info.printit, Sensitive=1
   Widget_Control, info.colorsID, Sensitive=1
ENDELSE

   ; Update the pixmap with histogram with no threshold lines.

XStretch_Histoplot, info.image, Background=info.backColor, Color=info.drawColor, $
   Max_Value=info.maxValue, _Extra=info.extra, WID=info.pixmap

      ; Put the info structure back into its storage location and then,
      ; out of here!

   Widget_Control, event.top, Set_UValue=info, /No_Copy
   RETURN
ENDIF ; thisEvent = UP


   ; Most of the action in this event handler occurs here while we are waiting
   ; for an UP event to occur. As long as we don't get it, keep erasing the
   ; old threshold line and drawing a new one.

   ; Get current window and scaling parameters in order.

WSet, info.histo_wid
!X.S = info.xs
!Y.S = info.ys
coord = Convert_Coord(event.x, event.y, /Device, /To_Data)

   ; Draw the "other" line on the pixmap (so you don't have to draw
   ; it all the time).

WSet, info.pixmap
CASE info.lineby OF
   'MIN': PlotS, [info.maxthresh, info.maxthresh],[info.ymin, info.ymax],  $
      Color=info.maxColor, Thick=3
   'MAX': PlotS, [info.minthresh, info.minthresh],[info.ymin, info.ymax],  $
      Color=info.minColor, Thick=3
ENDCASE

   ; Erase the old threshold line.

WSet, info.histo_wid
Device, Copy = [0, 0, info.pix_xsize, info.pix_ysize, 0, 0, info.pixmap]

   ; Draw the new line at the new coordinate. Make sure the coordinate
   ; is inside the plot and doesn't go over the other line.

CASE info.lineby OF
   'MIN': BEGIN
          coord(0) = coord(0) > (info.xmin + 1)
          coord(0) = coord(0) < (info.maxThresh - 1)
          END
   'MAX': BEGIN
          coord(0) = coord(0) > (info.minThresh + 1)
          coord(0) = coord(0) < (info.xmax - 1)
          END
ENDCASE

CASE info.lineby OF
   'MIN': PlotS, [coord(0), coord(0)],[info.ymin, info.ymax], $
       Color=info.minColor, Thick=3
   'MAX': PlotS, [coord(0), coord(0)],[info.ymin, info.ymax], $
       Color=info.maxColor, Thick=3
ENDCASE

   ; Put the info structure back into its storage location.

Widget_Control, event.top, Set_UValue=info, /No_Copy
END ; of XSTRETCH_MOVELINE **************************************************



PRO XSTRETCH_QUIT, event
Widget_Control, event.top, /Destroy
END ; of XSTRETCH_QUIT ******************************************************



PRO XSTRETCH_COLORS, event

Widget_Control, event.top, Get_UValue=info, /No_Copy

thisEvent = Tag_Names(event, /Structure_Name)
CASE thisEvent OF

   'WIDGET_BUTTON': BEGIN
       XColors, Group=event.top, NColors=info.ncolors, $
          NotifyID=[event.id, event.top]
       END
   'XCOLORS_LOAD': BEGIN
       Device, Get_Visual_Depth=thisDepth
       IF thisDepth GT 8 THEN BEGIN
          WSet, info.windex
          displayImage = BytScl(*info.image, Top=info.ncolors-1,  $
             Max=info.maxThresh, Min=info.minThresh)
          TVImage, displayImage, _Extra=info.extra
       ENDIF
       END
ENDCASE

Widget_Control, event.top, Set_UValue=info, /No_Copy

END ; of XSTRETCH_COLORS ****************************************************



PRO XSTRETCH_MAXVALUE, event

Widget_Control, event.top, Get_UValue=info, /No_Copy

   ; Get the new max value.

Widget_Control, event.id, Get_UValue=maxValue
info.maxValue = maxValue

   ; Update the histogram plot.

XStretch_Histoplot, info.image, Background=info.backColor, Color=info.drawColor, $
   Max_Value=info.maxValue, _Extra=info.extra, WID=info.histo_wid

   ; Draw threshold lines on the histogram plot.

WSet, info.histo_wid
PlotS, [info.minThresh, info.minThresh], [!Y.CRange(0), !Y.CRange(1)], $
   Color=info.minColor, Thick=3
PlotS, [info.maxThresh, info.maxThresh], [!Y.CRange(0), !Y.CRange(1)], $
   Color=info.maxColor, Thick=3

   ; Update the pixmap with histogram with no threshold lines.

XStretch_Histoplot, info.image, Background=info.backColor, Color=info.drawColor, $
   Max_Value=info.maxValue, _Extra=info.extra, WID=info.pixmap

Widget_Control, event.top, Set_UValue=info, /No_Copy

END ; of XSTRETCH_IMAGE_RESIZE **********************************************



PRO XSTRETCH_IMAGE_RESIZE, event

Widget_Control, event.top, Get_UValue=histoTLB
Widget_Control, histoTLB, Get_UValue=info, /No_Copy

Widget_Control, info.image_draw, Draw_XSize=event.x, Draw_YSize=event.y
WSet, info.windex
displayImage = BytScl(*info.image, Top=info.ncolors-1,  $
   Max=info.maxThresh, Min=info.minThresh)
TVImage, displayImage, _Extra=info.extra

Widget_Control, histoTLB, Set_UValue=info, /No_Copy

END ; of XSTRETCH_IMAGE_RESIZE **********************************************



PRO XSTRETCH_HISTOGRAM_RESIZE, event

Widget_Control, event.top, Get_UValue=info, /No_Copy

Widget_Control, info.histo_draw, Draw_XSize=event.x, Draw_YSize=event.y

   ; Draw the plot.

XStretch_Histoplot, info.image, Background=info.backColor, Color=info.drawColor, $
   Max_Value=info.maxValue, _Extra=info.extra, WID=info.histo_wid

   ; Put the same plot in the pixmap.

WDelete, info.pixmap
Window, /Free, XSize=event.x, YSize=event.y, /Pixmap
info.pixmap = !D.Window
info.pix_xsize = event.x
info.pix_ysize = event.y
Device, Copy=[0, 0, info.pix_xsize, info.pix_ysize, 0, 0, info.histo_wid]

   ; Save the scaling factors for calculating data coordinates.

info.xs = !X.S
info.ys = !Y.S

   ; Draw threshold lines on the histogram plot.

WSet, info.histo_wid
PlotS, [info.minThresh, info.minThresh], [!Y.CRange(0), !Y.CRange(1)], $
   Color=info.minColor, Thick=3
PlotS, [info.maxThresh, info.maxThresh], [!Y.CRange(0), !Y.CRange(1)], $
   Color=info.maxColor, Thick=3

Widget_Control, event.top, Set_UValue=info, /No_Copy

END ; of XSTRETCH_COLORS ****************************************************



PRO XSTRETCH_CLEANUP, tlb
Widget_Control, tlb, Get_UValue=info
IF N_Elements(info) NE 0 THEN BEGIN
   IF info.newPointer THEN Ptr_Free, info.image
   WDelete, info.pixmap
ENDIF
END ;---------------------------------------------------------------------



PRO XSTRETCH, theImage, Group_Leader=group, NColors=ncolors, $
   Max_Value=maxValue, Colortable=ctable, _EXTRA=extra

histXsize = 500
histYsize = 350

On_Error, 1
Device, Decomposed = 0

   ; Need an image?

IF N_Elements(theImage) EQ 0  THEN BEGIN
   file = Filepath(SubDir=['examples', 'data'], 'ctscan.dat')
   theImage = BytArr(256, 256)
   OpenR, lun, file, /GET_LUN
   ReadU, lun, theImage
   Free_LUN, lun
ENDIF

   ; Is image a pointer? If not, make it one.

IF Size(theImage, /TName) NE 'POINTER' THEN BEGIN
   image = Ptr_New(theImage)
   newPointer = 1
ENDIF ELSE BEGIN
   image = theImage
   newPointer = 0
ENDELSE

imgsize = Size(*image)
IF imgsize(0) NE 2 THEN $
   Message, 'First positional parameter must be a 2D image.'

xsize = imgsize(1)
ysize = imgsize(2)

  ; Default values for keywords.

IF N_Elements(maxValue) EQ 0 THEN BEGIN
   numPixels = N_Elements(*image)
   maxValue =  25000.0
   IF numPixels LE 65536L THEN maxValue = 5000.0
   IF numPixels GT 65536L AND numPixels LE 262144L THEN maxValue = 20000.0
ENDIF
IF N_Elements(extra) EQ 0 THEN extra = {TITLE:''}
IF N_Elements(ncolors) EQ 0 THEN BEGIN

   ; Check for availability of GIF files.

thisVersion = Float(!Version.Release)
IF thisVersion LT 5.4 THEN haveGif = 1 ELSE haveGIF = 0

      ; Find out how many colors you have.

   ncolors = !D.Table_Size - 4

   IF ncolors LT 24 THEN BEGIN
      Message, 'Not enough colors available to continue. Returning.'
      RETURN
   ENDIF

   minColor = ncolors
   maxColor = ncolors + 1
   backColor = ncolors + 2
   drawColor = ncolors + 3

ENDIF ELSE BEGIN

      ; We will scale to as many colors as we have, less 4 drawing colors.
      ; Must have at least 20 data colors.

   officialColors = !D.Table_Size
   ncolors = (ncolors-4) < (officialColors-4)
   IF ncolors LT 24 THEN BEGIN
      Message, 'Not enough colors available to continue. Returning.'
      RETURN
   ENDIF
   minColor = ncolors
   maxColor = ncolors + 1
   backColor = ncolors + 2
   drawColor = ncolors + 3

ENDELSE

   ; Create the histogram widget.

histo_tlb = Widget_Base(Row=1, Title='XStretch Histogram', $
   MBar=menubaseID, TLB_Size_Events=1, XOffset=50, YOffset=100)
histo_draw = Widget_Draw(histo_tlb, XSize=histXsize, YSize=histYsize, $
   Button_Events=1, Event_Pro='XStretch_Process_Events')
controlID = Widget_Button(menubaseID, Value='Controls', Event_Pro='XStretch_MaxValue')
saveAs = Widget_Button(controlID, Value='Save Image As...', Event_Pro="XStretch_SaveAs", /Menu)
dummy = Widget_Button(saveAs, Value='JPEG File', UValue='JPEG')
dummy = Widget_Button(saveAs, Value='TIFF File', UValue='TIFF')
IF havegif THEN dummy = Widget_Button(saveAs, Value='GIF File', UValue='GIF')
dummy = Widget_Button(saveAs, Value='PostScript File', UValue='PS')
printit = Widget_Button(controlID, Value='Print Image...', Event_Pro='XStretch_Print')
maxID = Widget_Button(controlID, Value='Max Pixel Density', /Menu, /Separator)
dummy = Widget_Button(maxID, Value='2000', UValue=2000.0)
dummy = Widget_Button(maxID, Value='5000', UValue=5000.0)
dummy = Widget_Button(maxID, Value='10000', UValue=10000.0)
dummy = Widget_Button(maxID, Value='20000', UValue=20000.0)
dummy = Widget_Button(maxID, Value='30000', UValue=30000.0)
dummy = Widget_Button(maxID, Value='50000', UValue=50000.0)
dummy = Widget_Button(maxID, Value='75000', UValue=75000.0)
colorsID = Widget_Button(controlID, Value='Image Colors...', $
   Event_Pro='XStretch_Colors')
quitter = Widget_Button(controlID, Value='Quit', $
   Event_Pro='XStretch_Quit', /Separator)
Widget_Control, histo_tlb, /Realize

   ; Create a pixmap window for moving and erasing the histogram
   ; threshold bars.

Window, Pixmap=1, XSize=histXsize, YSize=histYsize, /Free
pixmap = !D.Window

   ; Create an image window for displaying the image.

Widget_Control, histo_tlb, TLB_Get_Offset=offsets, TLB_Get_Size=sizes
xoff = offsets[0] + sizes[0] + 20
yoff = offsets[1]
aspect = Float(xsize)/ysize
IF xsize GT 512 OR ysize GT 512 THEN BEGIN
   IF xsize GT ysize THEN BEGIN
      xsize = 512
      ysize = Fix(xsize / aspect)
   ENDIF ELSE BEGIN
      ysize = 512
      xsize = ysize / aspect
   ENDELSE
ENDIF
image_tlb = Widget_Base(Row=1, Group=histo_tlb, Title='XStretch Image', $
   XOffSet=xoff, YOffSet=yoff, TLB_Size_Events=1)
image_draw = Widget_Draw(image_tlb, XSize=xsize, YSize=ysize, $
   Kill_Notify='XStretch_ImageWindowKilled', UValue=[saveAs, printit, colorsID])
Widget_Control, image_tlb, /Realize

   ; Get window index numbers for the draw widgets.

Widget_Control, image_draw, Get_Value=windex
Widget_Control, histo_draw, Get_Value=histo_wid

   ; Load a colortable if requested.

IF N_Elements(ctable) NE 0 THEN $
   LoadCt, 0 > ctable < 40, NColors=ncolors, /Silent ELSE $
   LoadCT, 0, NColors=ncolors, /Silent

   ; Load drawing colors.

TVLct, 255b, 255b, 0b, minColor    ; Yellow color.
TVLct, 0b, 255b, 0b, maxColor      ; Green color
TVLct, 70b, 70b, 70b, backColor    ; Charcoal color
TvLct, 255b, 255b, 255b, drawColor ; White color

   ; Get the current color table vectors for storage.

TVLCT, r, g, b, /Get

   ; Start with 2% linear stretch on both ends.

maxVal = Max(*image)
maxThresh = 0.98 * maxVal
minVal = Min(*image)
minThresh = minVal + (0.02 * maxVal)
XStretch_Histoplot, image, Background=backColor, Color=drawColor, $
   Max_Value=maxValue, _Extra=extra, WID=histo_wid

   ; Put the same plot in the pixmap.

WSet, pixmap
Device, Copy=[0, 0, histXsize, histYsize, 0, 0, histo_wid]

   ; Save the scaling factors for calculating data coordinates.

xs = !X.S
ys = !Y.S

WSet, histo_wid

   ; Draw threshold lines.

PlotS, [minThresh, minThresh], [!Y.CRange(0), !Y.CRange(1)], $
   Color=minColor, Thick=3
PlotS, [maxThresh, maxThresh], [!Y.CRange(0), !Y.CRange(1)], $
   Color=maxColor, Thick=3

   ; Display the image after thresholding.

WSet, windex
displayImage = BytScl(*image, Top=ncolors-1, Max=maxThresh, Min=minThresh)
TVImage, displayImage, _Extra=extra

   ; Calculate a value to tell you if you are "close" to a threshold line.

close = 0.05 * (maxval-minval)

  ; Make an info structure with all info to run the program.

info = {image:image, $           ; A pointer to the image data
        minThresh:minThresh, $   ; The minimum threshold
        maxThresh:maxThresh, $   ; The maximum threshold
        ncolors:ncolors, $       ; The number of colors
        minColor:minColor, $     ; The minimum drawing color index
        maxColor:maxColor, $     ; The maximum drawing color index
        backColor:backColor, $   ; The background drawing color index
        drawColor:drawColor, $   ; The plot drawing color index
        histo_wid:histo_wid, $   ; The histogram window index number
        histo_draw:histo_draw, $ ; The histogram draw widget ID.
        image_draw:image_draw, $ ; The image draw widget ID.
        maxValue:maxValue, $     ; The maximum value of the plot
        windex:windex, $         ; The image window index
        ymin:!Y.Crange(0), $     ; The ymin in data coordinates
        ymax:!Y.Crange(1), $     ; The ymax in data coordinates
        xmin:!X.Crange(0), $     ; The xmin in data coordinates
        xmax:!X.Crange(1), $     ; The xmax in data coordinates
        lineby:'MIN', $          ; The line you are close to.
        linex:minThresh, $       ; The x coordinate of line (data coords).
        pixmap:pixmap, $         ; The pixmap window index
        minval:minval, $         ; The minimum intensity value of the data
        maxval:maxval, $         ; The maximum intensity value of the data
        r:r, $                   ; Original red colors to restore.
        g:g, $                   ; Original green colors to restore.
        b:b, $                   ; Original blue colors to restore.
        extra:extra, $           ; The extra keywords for the Plot command.
        xs:xs, $                 ; Scaling x factors
        ys:ys, $                 ; Scaling y factors
        pix_xsize:histXsize, $   ; The X size of the pixmap.
        pix_ysize:histYsize, $   ; The Y size of the pixmap.
        newPointer:newPointer, $ ; A flag that indicates if we made a pointer or not.
        saveAs:saveAs, $         ; The SaveAs button widget identifier.
        printIt:printIt, $       ; The Print button widget identifier.
        colorsID:colorsID, $     ; The Image Colors button widget identifier.
        close:close}             ; A value to indicate closeness to line

   ; Save the info structure and bring the histogram window forward with SHOW.

Widget_Control, histo_tlb, Set_UValue=info, /No_Copy, /Show
Widget_Control, image_tlb, Set_UValue=histo_tlb
XManager, 'xstretch', histo_tlb, Group=group, /Just_Reg, /No_Block, $
   Event_Handler='XStretch_Histogram_Resize', Cleanup='XStretch_Cleanup'
XManager, 'xstretch_image', image_tlb, Event_Handler='XStretch_Image_Resize', /No_Block

END
