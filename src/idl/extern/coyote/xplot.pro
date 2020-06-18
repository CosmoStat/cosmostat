;+
; NAME:
;       XPLOT
;
; PURPOSE:
;       The purpose of this program is to demonstrate how to
;       create a line plot with axes and a title in the
;       new IDL 5 object graphics.
;
; AUTHOR:
;       FANNING SOFTWARE CONSULTING
;       David Fanning, Ph.D.
;       1645 Sheely Drive
;       Fort Collins, CO 80526 USA
;       Phone: 970-221-0438
;       E-mail: davidf@dfanning.com
;       Coyote's Guide to IDL Programming: http://www.dfanning.com
;
; CATEGORY:
;       Widgets, IDL 5 Object Graphics.
;
; CALLING SEQUENCE:
;       XPlot, x, y
;
; REQUIRED INPUTS:
;       x: A vector of input values used as the dependent data.
;
; OPTIONAL INPUTS
;       y: A vector of input values used as the dependent data.
;          If both x and y parameters are present, x is the independent data.
;
; OPTIONAL KEYWORD PARAMETERS:
;
;       COLORPRINT: Set this keyword to obtain printed output in color.
;       If not set (the default), printed output is in black and white.
;
;       EXACT: Set this keyword to a one- or two-element array to set exact axis
;       scaling for the axes. If Exact is a one-element array, both axes are
;       set to the same value. If Exact is a two-element array, the first
;       elements sets the X axis property and the second element sets the Y
;       axis property. For example, to set the X axis to exact scaling and
;       the Y axis to normal scaling, type:
;
;           IDL> x = Findgen(10)
;           IDL> XPlot, x, Sin(x), Exact=[1,0], XRange=[0,8.5]
;
;       _EXTRA: This keyword collects otherwise undefined keywords that are
;       passed to new Plot command. To some extent these are similar to the
;       old IDL Plot command. For example: Linestyle=2, Thick=3,
;       XRange=[-100,100], etc.
;
;       GROUP_LEADER: The group leader for this program. When the group leader
;       is destroyed, this program will be destroyed.
;
;       LANDSCAPE: Set this keyword if you are printing in landscape mode. The
;       default is Portrait mode. The Landscape keyword on the PRINTER object
;       is set, but not all printers will honor this keyword setting. If yours
;       does not, set Landscape mode in the Printer Setup dialog.
;
;       PSYM: The index of a plotting symbol to use on the plot. Integers 0-7
;       are valid values.
;
;       SYMSIZE: Sets the size of the symbols. By default, symbols are sized
;       so that they are 0.015 percent of the axis range.
;
;       VECTOR: Set this keyword if you want the printed output to be in
;       vector (as opposed to bitmap) form. This is faster, but not as accurate.
;
;       TITLE: A string used as the title of the plot.
;
;       XTITLE: A string used as the X title of the plot.
;
;       YTITLE: A string used as the Y title of the plot.
;
; COMMON BLOCKS:
;       None.
;
; SIDE EFFECTS:
;       None.
;
; RESTRICTIONS:
;       None.
;
; EXAMPLE:
;       To use this program, pass a 1D vector or vectors, like this:
;
;        IDL> XPlot, RandomU(seed, 11) * 9, YRange=[0, 10]
;
; MODIFICATION HISTORY:
;       Written by David Fanning, 13 June 97.
;       Modified axis font handling. 17 Sept 97. DWF.
;       Was not destroying all objects on exit. 12 Feb 98. DWF.
;       Changed IDLgrContainer to IDL_Container to fix 5.1 problems. 20 May 98. DWF.
;       Fixed a bug in the way symbols were (NOT!) sized. 11 May 99. DWF.
;       Added non-exact axis scaling. 12 May 99. DWF.
;       Fixed bug that changed data when calling with single parameter. 13 May DWF.
;       Added VECTOR, LANDSCAPE and COLORPRINT keywords and improved printing
;          capabilities. 16 Feb 2000. DWF.
;       Modified the EXACT keyword to accept values for X and Y axes
;          independently. 10 May 2000. DWF.
;       Updated for IDL 5.4. 13 June 2001. DWF.
;-
;
;###########################################################################
;
; LICENSE
;
; This software is OSI Certified Open Source Software.
; OSI Certified is a certification mark of the Open Source Initiative.
;
; Copyright © 1997-2000 Fanning Software Consulting.
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

    ; This is a utility routine to calculate the scaling vector
    ; required to position a vector of specified range at a
    ; specific position given in normalized coordinates. The
    ; scaling vector is given as a two-element array like this:
    ;
    ;   scalingVector = [translationFactor, scalingFactor]
    ;
    ; The scaling vector should be used with the [XYZ]COORD_CONV
    ; keywords of a graphics object or model. For example, if you
    ; wanted to scale an X axis into the data range of -0.5 to 0.5,
    ; you might type something like this:
    ;
    ;   xAxis->GetProperty, Range=xRange
    ;   xScale = Normalize(xRange, Position=[-0.5, 0.5])
    ;   xAxis, XCoord_Conv=xScale

IF (N_Elements(position) EQ 0) THEN position = [0.0, 1.0] ELSE $
    position=Float(position)
range = Float(range)

scale = [((position[0]*range[1])-(position[1]*range[0])) / $
    (range[1]-range[0]), (position[1]-position[0])/(range[1]-range[0])]

RETURN, scale
END
;-------------------------------------------------------------------------



FUNCTION XPlot_Aspect, aspectRatio, MARGIN=margin, WindowAspect=wAspectRatio

; This function calculates the correct aspect ratios for printing.

ON_ERROR, 2

   ; Check for aspect ratio parameter and possibilities.

IF N_PARAMS() EQ 0 THEN aspectRatio = 1.0

IF aspectRatio EQ 0 THEN BEGIN
   MESSAGE, 'Aspect Ratio of 0. Changing to 1...', /Informational
   aspectRatio = 1.0
ENDIF

s = SIZE(aspectRatio)
IF s(s(0)+1) NE 4 THEN $
   MESSAGE, 'Aspect Ratio is not a FLOAT. Take care...', /Informational

   ; Check for margins.

IF N_ELEMENTS(margin) EQ 0 THEN margin = 0.15

   ; Error checking.

IF margin LT 0 OR margin GE 0.5 THEN $
   MESSAGE, 'The MARGIN keyword value must be between 0.0 and 0.5.'

   ; Calculate the aspect ratio of the current window.

IF N_Elements(wAspectRatio) EQ 0 THEN wAspectRatio = FLOAT(!D.Y_VSIZE) / !D.X_VSIZE

   ; Calculate normalized positions in window.

IF (aspectRatio LE wAspectRatio) THEN BEGIN
   xstart = margin
   ystart = 0.5 - (0.5 - margin) * (aspectRatio / wAspectRatio)
   xend = 1.0 - margin
   yend = 0.5 + (0.5 - margin) * (aspectRatio / wAspectRatio)
ENDIF ELSE BEGIN
   xstart = 0.5 - (0.5 - margin) * (wAspectRatio / aspectRatio)
   ystart = margin
   xend = 0.5 + (0.5 - margin) * (wAspectRatio / aspectRatio)
   yend = 1.0 - margin
ENDELSE

position = [xstart, ystart, xend, yend]

RETURN, position
END
;-------------------------------------------------------------------------



PRO XPlot_Output, event

   ; This event handler creates GIF and JPEG files.

Widget_Control, event.top, Get_UValue=info, /No_Copy

   ; Get a snapshop of window contents. (TVRD equivalent.)

info.thisWindow->GetProperty, Image_Data=snapshot

   ; What kind of file is wanted?

Widget_Control, event.id, GET_UValue=whichFileType
CASE whichFileType OF

   'GIF': BEGIN

         ; Because we are using RGB color for this model, we have
         ; a 3-m-n array. Use Color_Quan to create a 2D image and
         ; appropriate color tables for the GIF file.

      image2D = Color_Quan(snapshot, 1, r, g, b)
      filename = Dialog_Pickfile(/Write, File='xplot.gif')
      IF filename NE '' THEN Write_GIF, filename, image2d, r, g, b
      END

   'JPEG': BEGIN

      filename = Dialog_Pickfile(/Write, File='xplot.jpg')
      IF filename NE '' THEN Write_JPEG, filename, snapshot, True=1
      END


   'TIFF': BEGIN

      filename = Dialog_Pickfile(/Write, File='xplot.tif')
      IF filename NE '' THEN BEGIN

         ; TIFF files should have their Y direction reversed for
         ; compatibility with most other software.

         Write_TIFF, filename, Reverse(snapshot,3)
      ENDIF
      END

ENDCASE

    ;Put the info structure back.

Widget_Control, event.top, Set_UValue=info, /No_Copy
END
;-------------------------------------------------------------------



PRO XPlot_Exit, event

   ; Exit the program.

Widget_Control, event.top, /Destroy
END
;-------------------------------------------------------------------



PRO XPlot_Printing, event

   ; Printer output handled here.

Widget_Control, event.top, Get_UValue=info, /No_Copy

   ; Does the user really want to print?

print_it = Dialog_PrinterSetup(info.thisPrinter)
IF NOT print_it THEN BEGIN
   Widget_Control, event.top, Set_UValue=info, /No_Copy
   RETURN
ENDIF

   ; Find out the current colors of all the objects.

info.plotView->GetProperty, Color=backgroundColor
info.thisPlot->GetProperty, Color=plotColor
info.xAxis1->GetProperty, Color=axisColor
info.thisSymbol->GetProperty, Color=symbolColor

   ; Change colors to black and white for printing.

IF NOT info.colorprint THEN BEGIN
   info.plotView->SetProperty, Color=[255,255,255]
   info.thisPlot->SetProperty, Color=[0,0,0]
   info.xAxis1->SetProperty, Color=[0,0,0]
   info.yAxis1->SetProperty, Color=[0,0,0]
   info.xAxis2->SetProperty, Color=[0,0,0]
   info.yAxis2->SetProperty, Color=[0,0,0]
   info.thisSymbol->SetProperty, Color=[0,0,0]
   info.plotTitle->SetProperty, Color=[0,0,0]
ENDIF

   ; I want the output on the page to have the same aspect ratio
   ; (ratio of height to width) as I see in the display window.

info.thisWindow->GetProperty, Dimensions=wdims
info.thisPrinter->GetProperty, Dimensions=pdims
plotAspect = Float(wdims[1]) / wdims[0]
windowAspect = Float(pdims[1]) / pdims[0]
position = XPlot_Aspect(plotAspect, WindowAspect=windowAspect, Margin=0)
info.plotView->SetProperty, Dimensions=[position[2]-position[0], position[3]-position[1]], $
   Location=[position[0], position[1]], Units=3

   ; Print it. May take a little time. Alert the user.

Widget_Control, Hourglass=1
info.thisPrinter->Draw, info.plotView, Vector=info.vector
info.thisPrinter->NewDocument
Widget_Control, Hourglass=0

   ; Put everything back the way it was.

IF NOT info.colorprint THEN BEGIN
   info.plotView->SetProperty, Color=backgroundColor
   info.thisPlot->SetProperty, Color=plotColor
   info.xAxis1->SetProperty, Color=axisColor
   info.yAxis1->SetProperty, Color=axisColor
   info.xAxis2->SetProperty, Color=axisColor
   info.yAxis2->SetProperty, Color=axisColor
   info.thisSymbol->SetProperty, Color=symbolColor
   info.plotTitle->SetProperty, Color=axisColor
ENDIF

info.PlotView->SetProperty, Location=[0,0], Dimensions=[0,0]

   ; Put the info structure back.

Widget_Control, event.top, Set_UValue=info, /No_Copy
END
;-------------------------------------------------------------------



PRO Xplot_Cleanup, id

    ; Come here when the widget dies. Free all the program
    ; objects, pointers, pixmaps, etc. and release memory.

Widget_Control, id, Get_UValue=info
IF N_Elements(info) NE 0 THEN Obj_Destroy, info.thisContainer
END
;---------------------------------------------------------------------



PRO XPlot_Draw_Widget_Events, event

    ; This event handler handles draw widget expose events.

Widget_Control, event.top, Get_UValue=info, /No_Copy

    ; Draw the graphic.

info.thisWindow->Draw, info.plotView

    ;Put the info structure back.

Widget_Control, event.top, Set_UValue=info, /No_Copy
END
;------------------------------------------------------------------------



PRO XPlot_Resize_Events, event

    ; This event handler handles TLB resize events.

Widget_Control, event.top, Get_UValue=info, /No_Copy

    ; Resize the draw widget.

info.thisWindow->SetProperty, Dimension=[event.x, event.y]

    ; Redisplay the graphic.

info.thisWindow->Draw, info.plotView

    ;Put the info structure back.

Widget_Control, event.top, Set_UValue=info, /No_Copy
END
;------------------------------------------------------------------------



PRO XPlot_Linestyle, event

    ; This event handler handles linesytle change events.

Widget_Control, event.top, Get_UValue=info, /No_Copy

    ; Get the requested linestyle.

Widget_Control, event.id, Get_UValue=thisLineStyle

    ; Change it.

info.thisPlot->SetProperty, LineStyle=thisLineStyle

    ; Redisplay the graphic.

info.thisWindow->Draw, info.plotView

    ;Put the info structure back.

Widget_Control, event.top, Set_UValue=info, /No_Copy
END
;------------------------------------------------------------------------



PRO XPlot_Symbol, event

    ; This event handler handles symbol change events.

Widget_Control, event.top, Get_UValue=info, /No_Copy

    ; Get the requested symbol.

Widget_Control, event.id, Get_UValue=thisSymbol

    ; Change it.

info.thisPlot->GetProperty, Symbol=symbolObject
symbolObject->SetProperty, Data=thisSymbol

    ; Redisplay the graphic.

info.thisWindow->Draw, info.plotView

    ;Put the info structure back.

Widget_Control, event.top, Set_UValue=info, /No_Copy
END
;------------------------------------------------------------------------



PRO XPlot_SymbolSize, event

    ; This event handler handles symbol size change events.

Widget_Control, event.top, Get_UValue=info, /No_Copy

    ; Get the requested symbol size.

Widget_Control, event.id, Get_UValue=thisSize

    ; Change it.

info.thisPlot->GetProperty, Symbol=symbolObject
symbolObject->SetProperty, Size=info.symbolSize * thisSize

    ; Redisplay the graphic.

info.thisWindow->Draw, info.plotView

    ;Put the info structure back.

Widget_Control, event.top, Set_UValue=info, /No_Copy
END
;------------------------------------------------------------------------


PRO XPlot, xx, yy, _Extra=extra, PSym=psym, Title=title, SymSize=symSize, $
   Group_Leader=group, XTitle=xtitle, YTitle=ytitle, Exact=exact, $
   ColorPrint=colorprint, Vector=vector, Landscape=landscape

   ; New printer functionality requires IDL 5.3 or higher.

IF Float(!Version.Release) LT 5.3 THEN BEGIN
   ok = Dialog_Message('Program functionality requires IDL 5.3 or higher. Returning...')
   RETURN
ENDIF

    ; Check to be sure at least one parameter is present.

np =  N_Params()
CASE np OF
    0: BEGIN
       Print, 'Using fake data in XPLOT...'
       y = FIndGen(101)
       y = Sin(y/5) / Exp(y/50)
       x = IndGen(N_Elements(y))
       xtitle = 'Time'
       ytitle = 'Signal Stength'
       title = 'Example Line Plot'
       END
    1: BEGIN
       y = xx
       x = IndGen(N_Elements(y))
       END
    ELSE:
ENDCASE

   ; Make sure no data parameters change.

IF N_Elements(x) EQ 0 THEN x = xx
IF N_Elements(y) EQ 0 THEN y = yy

    ; Check keyword parameters.

IF N_Elements(psym) EQ 0 THEN psym = 0
IF N_Elements(title) EQ 0 THEN title = ''
IF N_Elements(symsize) EQ 0 THEN symsize = 1.0
IF N_Elements(xtitle) EQ 0 THEN xtitle = ''
IF N_Elements(ytitle) EQ 0 THEN ytitle = ''
colorprint = Keyword_Set(colorprint)
vector = Keyword_Set(vector)
landscape = Keyword_Set(landscape)
CASE N_Elements(exact) OF
   0: exact = [0,0]
   1: exact = Replicate(exact, 2)
   2:
   ELSE: BEGIN
      ok = Dialog_Message('Exact keyword contains too many elements. Returning...')
      RETURN
      ENDCASE
ENDCASE

   ; Create title objects for the axes. Color them yellow.

xTitle = Obj_New('IDLgrText', xtitle, Color=[255,255,0])
yTitle = Obj_New('IDLgrText', ytitle, Color=[255,255,0])

    ; Make a symbol object. Color symbols cyan.

thisSymbol = Obj_New('IDLgrSymbol', psym, Color=[0, 255, 255])

    ; Make a font object.

helvetica10pt = Obj_New('IDLgrFont', 'Helvetica', Size=10)

    ; Create a plot object. The plot will be in the coordinate
    ; space 0->1. The view will be in the range -0.35->1.25 so
    ; that the plot axis annotation will be visable. Make the plot
    ; a green color.

thisPlot = Obj_New("IDLgrPLOT", x, y, _Extra=extra, $
   Color=[0,255,0], Symbol=thisSymbol, Thick=2)

    ; Get the data ranges from the Plot Object.

thisPlot->GetProperty, XRange=xrange, YRange=yrange

    ; Create plot box style axes. Make the axes yellow.
    ; The large values in the LOCATION keyword indicates which
    ; values are NOT used. The axes text is set to Helvetica
    ; 10 point font.

xAxis1 = Obj_New("IDLgrAxis", 0, Color=[255,255,0], Ticklen=0.025, $
    Minor=4, Range=xrange, Title=xtitle, $
    Location=[1000, 0 ,0], Exact=exact[0])
xAxis1->GetProperty, Ticktext=xAxisText
xAxisText->SetProperty, Font=helvetica10pt

xAxis2 = Obj_New("IDLgrAxis", 0, Color=[255,255,0], Ticklen=0.025, $
    Minor=4, /NoText, Range=xrange, TickDir=1, $
    Location=[1000, 1, 0], Exact=exact[0])

yAxis1 = Obj_New("IDLgrAxis", 1, Color=[255,255,0], Ticklen=0.025, $
    Minor=4, Title=ytitle, Range=yrange, $
    Location=[0, 1000, 0], Exact=exact[1])
yAxis1->GetProperty, Ticktext=yAxisText
yAxisText->SetProperty, Font=helvetica10pt

yAxis2 = Obj_New("IDLgrAxis", 1, Color=[255,255,0], Ticklen=0.025, $
    Minor=4, /NoText, Range=yrange, TickDir=1, $
    Location=[1, 1000, 0], Exact=exact[1])

    ; Because we may not be using exact axis ranging, the axes
    ; may extend further than the xrange and yrange. Get the
    ; actual axis range so that the plot, etc. can be scaled
    ; appropriately.

xAxis1->GetProperty, CRange=xrange
yAxis1->GetProperty, CRange=yrange

    ; Set up the scaling so that the axes for the plot and the
    ; plot data extends from 0->1 in the X and Y directions.

xs = Normalize(xrange)
ys = Normalize(yrange)

    ; Scale the plot data and axes into 0->1.

thisPlot->SetProperty, XCoord_Conv=xs, YCoord_Conv=ys
xAxis1->SetProperty, XCoord_Conv=xs
xAxis2->SetProperty, XCoord_Conv=xs
yAxis1->SetProperty, YCoord_Conv=ys
yAxis2->SetProperty, YCoord_Conv=ys

    ; Size the symbols appropriately for the plot.

xSymSize = (xrange[1] - xrange[0]) * 0.015 * symSize
ySymSize = (yrange[1] - yrange[0]) * 0.015 * symSize
IF Obj_Valid(thisSymbol) THEN thisSymbol->SetProperty, Size=[xSymSize, ySymSize]

  ; Create a plot title. Center it at a location above the plot.

helvetica14pt = Obj_New('IDLgrFont', 'Helvetica', Size=14)
plotTitle = Obj_New('IDLgrText', title, Color=[255,255,0], $
   Location=[0.5, 1.05, 0.0], Alignment=0.5, Font=helvetica14pt)

    ; Create a plot model and add axes, plot, and plot title to it.

plotModel = Obj_New('IDLgrModel')
plotModel->Add, thisPlot
plotModel->Add, xAxis1
plotModel->Add, xAxis2
plotModel->Add, yAxis1
plotModel->Add, yAxis2
plotModel->Add, plotTitle

    ; Create a view and add the plot model to it. Notice that the view
    ; is larger than the 0->1 plot area to accomodate axis annotation.
    ; The view will have a gray background.

plotView = Obj_New('IDLgrView', Viewplane_Rect=[-.35, -.35, 1.6, 1.6], $
   Location=[0,0], Color=[80,80,80])
plotView->Add, plotModel

   ; Check for availability of GIF files.

thisVersion = Float(!Version.Release)
IF thisVersion LT 5.4 THEN haveGif = 1 ELSE haveGIF = 0

    ; Create the widgets for this program.

tlb = Widget_Base(Column=1, Title='Resizeable Line Plot Example', $
   TLB_Size_Events=1, MBar=menubase)

    ; Create FILE menu buttons for printing and exiting.

filer = Widget_Button(menubase, Value='File', /Menu)
b = Widget_Button(filer, Value='Print', $
   Event_Pro='XPlot_Printing', UValue='PRINT')
b = Widget_Button(filer, /Separator, Value='Exit', $
   Event_Pro='XPlot_Exit')

    ; Create OUTPUT menu buttons for formatted output files.

output = Widget_Button(menubase, Value='Output')
IF havegif THEN b = Widget_Button(output, Value='GIF File', $
   UValue='GIF', Event_Pro='XPlot_Output')
b = Widget_Button(output, Value='JPEG File', $
   UValue='JPEG', Event_Pro='XPlot_Output')
b = Widget_Button(output, Value='TIFF File', $
   UValue='TIFF', Event_Pro='XPlot_Output')

    ; Create PROPERTIES menu buttons for plot properties.

propertiesID = Widget_Button(menubase, Value='Properties')

linestyleID = Widget_Button(propertiesID, Value='Line Style', $
   Event_Pro='XPlot_Linestyle', /Menu)
b = Widget_Button(linestyleID, Value='Solid', UValue=0)
b = Widget_Button(linestyleID, Value='Dot', UValue=1)
b = Widget_Button(linestyleID, Value='Dash', UValue=2)
b = Widget_Button(linestyleID, Value='Dash Dot', UValue=3)
b = Widget_Button(linestyleID, Value='Dash Dot Dot Dot', UValue=4)
b = Widget_Button(linestyleID, Value='Long Dash', UValue=5)
b = Widget_Button(linestyleID, Value='No Line', UValue=6)

symbolID = Widget_Button(propertiesID, Value='Symbol', $
   Event_Pro='XPlot_Symbol', /Menu)
b = Widget_Button(symbolID, Value='No Symbol', UValue=0)
b = Widget_Button(symbolID, Value='Plus Sign', UValue=1)
b = Widget_Button(symbolID, Value='Asterisk', UValue=2)
b = Widget_Button(symbolID, Value='Period', UValue=3)
b = Widget_Button(symbolID, Value='Diamond', UValue=4)
b = Widget_Button(symbolID, Value='Triangle', UValue=5)
b = Widget_Button(symbolID, Value='Square', UValue=6)
b = Widget_Button(symbolID, Value='X', UValue=7)

symbolSizeID = Widget_Button(propertiesID, Value='Symbol Size', $
   Event_Pro='XPlot_SymbolSize', /Menu)
b = Widget_Button(symbolSizeID, Value='0.25', UValue=0.25)
b = Widget_Button(symbolSizeID, Value='0.50', UValue=0.50)
b = Widget_Button(symbolSizeID, Value='0.75', UValue=0.75)
b = Widget_Button(symbolSizeID, Value='1.00', UValue=1.00)
b = Widget_Button(symbolSizeID, Value='1.25', UValue=1.25)
b = Widget_Button(symbolSizeID, Value='1.50', UValue=1.50)
b = Widget_Button(symbolSizeID, Value='1.75', UValue=1.75)
b = Widget_Button(symbolSizeID, Value='2.00', UValue=2.00)

    ; Create the draw widget. Use RGB color model. Be sure to set
    ; the Expose_Events keyword so that graphics are redisplayed
    ; properly. (This implies that RETAIN=0. The Graphics_Level
    ; keyword makes it a window object.

drawID = Widget_Draw(tlb, XSize=400, YSize=400, Color_Model=0, $
   Graphics_Level=2, Expose_Events=1, Retain=0, $
   Event_Pro='XPlot_Draw_Widget_Events')

    ; Realize the widgets and get the window object.

Widget_Control, tlb, /Realize
Widget_Control, drawID, Get_Value=thisWindow

    ; Display the plot in the window.

thisWindow->Draw, plotView

   ; Get a printer object for this graphic.

thisPrinter = Obj_New('IDLgrPrinter', Landscape=landscape)

   ; Create a container object to hold all the other
   ; objects. This will make it easy to free all the
   ; objects when we are finished with the program.

thisContainer = Obj_New('IDL_Container')
thisContainer->Add, thisWindow
thisContainer->Add, plotView
thisContainer->Add, thisPrinter
thisContainer->Add, helvetica10pt
thisContainer->Add, helvetica14pt
thisContainer->Add, xaxis1
thisContainer->Add, xaxis2
thisContainer->Add, yaxis1
thisContainer->Add, yaxis2
thisContainer->Add, plotTitle
thisContainer->Add, xTitle
thisContainer->Add, yTitle
thisContainer->Add, thisSymbol

   ; Create an info structure to hold program information.

info = { thisContainer:thisContainer, $     ; The container object.
         thisPlot:thisPlot, $               ; The plot object.
         symbolSize:[xSymSize, ySymSize], $ ; The default symbol size.
         thisWindow:thisWindow, $           ; The window object.
         thisSymbol:thisSymbol, $           ; The symbol object.
         plotTitle:plotTitle, $             ; The plot title object.
         xTitle:xTitle, $                   ; The X axis title object.
         yTitle:yTitle, $                   ; The Y axis title object.
         plotView:plotView, $               ; The view that will be rendered.
         xaxis1:xaxis1, $                   ; The X axis object.
         xaxis2:xaxis2, $                   ; The X axis object.
         yaxis1:yaxis1, $                   ; The Y axis object.
         yaxis2:yaxis2, $                   ; The Y axis object.
         colorprint:colorprint, $           ; A flag for color printing.
         vector:vector, $                   ; A flag for vector printing.
         landscape:landscape, $             ; A flag for landscape printing.
         thisPrinter:thisPrinter }          ; The printer object.

    ; Put the info stucture in the TLB.

Widget_Control, tlb, Set_UValue=info, /No_Copy

    ; Start the event loop.

XManager, 'xplot', tlb, Cleanup='XPlot_Cleanup', Group_Leader=group, $
   Event_Handler='XPlot_Resize_Events', /No_Block

END
;------------------------------------------------------------------------


