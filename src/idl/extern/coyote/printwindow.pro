;+
; NAME
;   PRINTWINDOW
;
;    This program sends the contents of the specified
;    window to the default printer. The current window
;    is used if a window index number is not provided.
;
; POSITIONAL PARAMETERS:
;
;   WID       The window index number of the window to send to the
;             printer. !D.Window used by default.
;
; KEYWORD PARAMETERS:
;
;   LANDSCAPE  If this keyword is set, the output is in Landscape
;              mode. Otherwise, Portrait mode is used.
;
;   PAGESIZE: Set this keyword to a string indicating the type
;             of PostScript page size you want. Current values are "LETTER",
;             "LEGAL", and "A4". Default is "LETTER".
;
; Written by David Fanning based on previous PRINT_IT program. 29 July 2000.
;-

FUNCTION PrintWindow_Error, theMessage, Traceback=traceback, NoName=noName

On_Error, 2

   ; Check for presence and type of message.

IF N_Elements(theMessage) EQ 0 THEN theMessage = !Error_State.Msg
s = Size(theMessage)
messageType = s[s[0]+1]
IF messageType NE 7 THEN BEGIN
   Message, "The message parameter must be a string.", _Extra=extra
ENDIF

   ; Get the call stack and the calling routine's name.

Help, Calls=callStack
callingRoutine = (Str_Sep(StrCompress(callStack[1])," "))[0]

   ; Are widgets supported? Doesn't matter in IDL 5.3 and higher.

widgetsSupported = ((!D.Flags AND 65536L) NE 0) OR Float(!Version.Release) GE 5.3
IF widgetsSupported THEN BEGIN
   IF Keyword_Set(noName) THEN answer = Dialog_Message(theMessage, _Extra=extra) ELSE BEGIN
      IF StrUpCase(callingRoutine) EQ "$MAIN$" THEN answer = Dialog_Message(theMessage, _Extra=extra) ELSE $
         answer = Dialog_Message(StrUpCase(callingRoutine) + ": " + theMessage)
   ENDELSE
ENDIF ELSE BEGIN
      Message, theMessage, /Continue, /NoPrint, /NoName, /NoPrefix, _Extra=extra
      Print, '%' + callingRoutine + ': ' + theMessage
      answer = 'OK'
ENDELSE

   ; Provide traceback information if requested.

IF Keyword_Set(traceback) THEN BEGIN
   Help, /Last_Message, Output=traceback
   Print,''
   Print, 'Traceback Report from Error_Message:'
   Print, ''
   FOR j=0,N_Elements(traceback)-1 DO Print, "     " + traceback[j]
ENDIF

RETURN, answer
END ;-----------------------------------------------------------------------------------


PRO PrintWindow, wid, Landscape=landscape, PageSize=pageSize

   ; Check parameters.

IF N_Params() EQ 0 THEN wid = !D.Window
landscape = Keyword_Set(landscape)
IF N_Elements(pagesize) EQ 0 THEN pagesize = 'LETTER' $
   ELSE pagesize = StrUpCase(pagesize)
CASE pagesize OF
   'LETTER': BEGIN
      shortside = 8.5
      longside = 11.0
      ENDCASE
   'LEGAL': BEGIN
      shortside = 8.5
      longside = 14.0
      ENDCASE
    'A4': BEGIN
      shortside = 8.27
      longside = 11.7
      ENDCASE
    ELSE: BEGIN
      Message, 'Unknown page size. Using LETTER...', /Informational
      shortside = 8.5
      longside = 11.0
      ENDCASE
ENDCASE

   ; Are we running on a 24-bit device?

Catch, theError
IF theError NE 0 THEN BEGIN
   theDepth = 8
   GOTO, testDepth
ENDIF
Device, Get_Visual_Depth=theDepth
testDepth:
IF theDepth GT 8 THEN truecolor = 1 ELSE truecolor = 0

   ; Main error handler.

Catch, theError
IF theError NE 0 THEN BEGIN
   Catch, /Cancel
   ok = PrintWindow_Error(!Error_State.Msg)
   RETURN
ENDIF

   ; Valid window ID?

IF wid LT 0 THEN Message, 'No graphics window currently open.', /NoName

   ; Get information about printer.

ok = Dialog_PrinterSetup()
IF NOT ok THEN RETURN

   ; Make the window current. Get contents.

thisWindow = !D.Window
WSet, wid
contents = TVRD(True=truecolor)

   ; Need a window on printer with same aspect ratio as current window.
   ; Your fudge factor (to account for offsets calculated from the printable
   ; edge of the page rather than the real edge of the page) may be different.

keywords = PSWindow(/Printer, Landscape=landscape, Fudge=0.25, PageSize=pagesize)

   ; Change the current device to PRINTER. Copy color table.

thisDevice = !D.Name

Set_Plot, 'PRINTER', /Copy

   ; Reset the PRINTER for proper calculations. Landscape mode
   ; has to be set separately from sizes. I don't know why.

Device, Scale_Factor=1, Portrait=1
Device, Landscape = keywords.landscape

   ; Configure the printer.

Device, _Extra=keywords

   ; Print it.

TV, contents, True=truecolor, XSize=keywords.xsize, $
   YSize=keywords.ysize, Inches=keywords.inches
Plots, [0,0,1,1,0], [0,1,1,0,0], /normal
Device, /Close_Document

   ; Clean up.

Set_Plot, thisDevice
WSet, thisWindow

END
