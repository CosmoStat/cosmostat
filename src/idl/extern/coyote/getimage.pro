;+
; NAME:
;       GETIMAGE
;
; PURPOSE:
;       The purpose of this function is to allow the user to open either
;       regular or XDR binary image files of two or three dimensions.
;
; CATEGORY:
;       Widgets, File I/O.
;
; CALLING SEQUENCE:
;       image = GETIMAGE(filename)
;
; OPTIONAL INPUTS:
;       filename: The name of the file to open for reading.
;
; OPTIONAL KEYWORD PARAMETERS:
;
;       CANCEL: An output variable that can be set to a named variable.
;       The value of the return variable will be 1 if the user clicked
;       the "Cancel" button or if there was a problem reading the file.
;
;       CATCH: Set this keyword to 0 if you wish to turn error catching OFF.
;
;       DIRECTORY: The name of the directory the file is located in. By
;       default the program looks in the "coyote" directory under the
;       main IDL directory, if one exists. Otherwise, it defaults to the
;       current directory.
;
;       FRAMES: The 3rd dimension of a 3D data set. Defaults to 0.
;
;       HEADER: The size of any header information in the file in BYTES.
;       Default is 0.
;
;       PARENT: The group leader for this widget program. The PARENT is
;       required if GETIMAGE is called from another widget program in order
;       to make this program a MODAL widget program.
;
;       XDR: Set this keyword if the binary file is of XDR type.
;
;       XOFFSET: This is the X offset of the program on the display. The
;       program will be placed approximately in the middle of the display
;       by default.
;
;       XSIZE: The size of the 1st dimension of the data.
;
;       YOFFSET: This is the Y offset of the program on the display. The
;       program will be placed approximately in the middle of the display
;       by default.
;
;       YSIZE: The size of the 2nd dimension of the data.
;
; COMMON BLOCKS:
;       None.
;
; SIDE EFFECTS:
;       A "CANCEL" operation is indicated by a 0 return value.
;       Any error in reading the file results in a 0 return value.
;
; RESTRICTIONS:
;       None.
;
; EXAMPLE:
;       To load the image "galaxy.dat" in the $IDL/examples/data
;       directory, type:
;
;       image = GETIMAGE('galaxy.dat', DIRECTORY=!DIR + '/examples/data', $
;          XSIZE=256, YSIZE=256, Cancel=cancelled, Parent=event.top)
;       IF NOT cancelled THEN TV, image
;
; MODIFICATION HISTORY:
;       Written by: David Fanning, 3 February 96.
;       Fixed bug that prevented reading INTEGER data. 19 Dec 96.
;       Modifed program for IDL 5 MODAL operation. 19 Oct 97.
;       Added CANCEL keyword. 27 Oct 97. DWF.
;       Fixed CANCLE keyword spelling. Sigh... 29 JUN 98. DWF.
;       Added COYOTE_FIELD, improved appearance. 19 NOV 99. DWF.
;       Updated with latest version of COYOTE_FIELD. 18 FEB 2000. DWF.
;       Added CATCH keyword so the program will break when I want
;       it to. :-) 18 MAR 2000. DWF.
;       Added GROUP_LEADER keyword, which is synonymous with PARENT. 31 MAR 2000. DWF.
;       Updated obsolete PICKFILE call to DIALOG_PICKFILE. 17 JAN 2001. DWF.
;-
;
;###########################################################################
;
; LICENSE
;
; This software is OSI Certified Open Source Software.
; OSI Certified is a certification mark of the Open Source Initiative.
;
; Copyright © 2000 Fanning Software Consulting.
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
;
; NAME:
;   COYOTE_FIELD
;
; PURPOSE:
;
;   The purpose of this compound widget is to provide an alternative
;   to the CW_FIELD widget offered in the IDL distribution. What has
;   always annoyed me about CW_FIELD is that the text widgets do not
;   look editable to the users on Windows platforms. This program
;   corrects that deficiency and adds some features that I think
;   would be helpful. For example, you can now assign an event handler
;   to the compound widget.
;
; AUTHOR:
;   FANNING SOFTWARE CONSULTING
;   David Fanning, Ph.D.
;   2642 Bradbury Court
;   Fort Collins, CO 80521 USA
;   Phone: 970-221-0438
;   E-mail: davidf@dfanning.com
;   Coyote's Guide to IDL Programming: http://www.dfanning.com/
;
; CATEGORY:
;
;   General programming.
;
; CALLING SEQUENCE:
;
;   fieldID = COYOTE_Field(parent, Title='X Size: ", Value=256, /IntegerValue)
;
; INPUT PARAMETERS:
;
;   parent -- The parent widget ID of the compound widget. Required.
;
; INPUT KEYWORDS:
;
;   Column -- Set this keyword to have the Label Widget above the Text Widget.
;   CR_Only -- Set this keyword if you only want Carriage Return events. If
;              this keyword is not set, all events are returned. No events
;              are returned unless the EVENT_PRO or EVENT_FUNC keywords are used.
;   DoubleValue -- Set this keyword if you want DOUBLE values returned.
;   Decimal -- Set this keyword to the number of digits to the right of the decimal
;              point in FLOATVALUE and DOUBLEVALUE numbers.
;   Digits -- Set this keyword to the number of digits permitted in INTERGERVALUE and LONGVALUE numbers.
;   Event_Func -- Set this keyword to the name of an Event Function. If this
;                 keyword is undefined and the Event_Pro keyword is undefined,
;                 all compound widget events are handled internally and not
;                 passed on to the parent widget.
;   Event_Pro -- Set this keyword to the name of an Event Procedure. If this
;                keyword is undefined and the Event_Func keyword is undefined,
;                all compound widget events are handled internally and not
;                passed on to the parent widget.
;   FieldFont -- The font name for the text in the Text Widget.
;   FloatValue -- Set this keyword for FLOAT values.
;   Frame -- Set this keyword to put a frame around the compound widget.
;   IntegerValue -- Set this keyword for INTEGER values.
;   LabelFont -- The font name for the text in the Label Widget.
;   LabelSize -- The X screen size of the Label Widget.
;   LongValue -- Set this keyword for LONG values.
;   Row=row -- Set this keyword to have the Label beside the Text Widget. (The default.)
;   Scr_XSize -- The X screen size of the compound widget.
;   Scr_YSize -- The Y screen size of the compound widget.
;   StringValue -- Set this keyword for STRING values. (The default.)
;   Title -- The text to go on the Label Widget.
;   UValue -- A user value for any purpose.
;   Value -- The "value" of the compound widget.
;   XSize -- The X size of the Text Widget.
;
; COMMON BLOCKS:
;
;   None.
;
; RESTRICTIONS:
;
;   None.
;
; EVENT STRUCTURE:
;
;   All events are handled internally unless either the Event_Pro or Event_Func
;   keywords are used to assign an event handler to the compound widget. By
;   default all events generated by the text widget are passed to the assigned
;   event handler. If you wish to receive only Carriage Return events, set the
;   CR_Only keyword.
;
;   event = { COYOTE_FIELD, $      ; The name of the event structure.
;             ID: 0L, $            ; The ID of the compound widget's top-level base.
;             TOP: 0L, $           ; The widget ID of the top-level base of the hierarchy.
;             HANDLER: 0L, $       ; The event handler ID. Filled out by IDL.
;             Value: Ptr_New(), $  ; A pointer to the widget value.
;             Type:""              ; A string indicating the type of data in the VALUE field.
;           }                      ; Values are "INT", "LONG", "FLOAT", "DOUBLE", or "STRING".
;
; EXAMPLE:
;
;   An example program is provided at the end of the COYOTE_FIELD code. To run it,
;   type these commands:
;
;      IDL> .Compile COYOTE_Field
;      IDL> Example
;
; MODIFICATION HISTORY:
;
;   Written by: David Fanning, 17 NOV 1999.
;   Added check to make float and double values finite. 18 NOV 1999. DWF.
;   Fixed a bug when selecting and deleting all numerical text. 19 NOV 1999. DWF.
;   Added DECIMAL and DIGITS keywords. 2 Jan 2000. DWF.
;   Added the POSITIVE keyword. 12 Jan 2000. DWF.
;   Fixed a few minor bugs with delete and digits. 12 Jan 2000. DWF.
;   Made GET_VALUE function return pointer to data, instead of data. 12 Jan 2000. DWF.
;



Function COYOTE_Field_ReturnValue, inputValue, dataType

; This utility routine takes a string and turns it into a number,
; depending upon the required data type.

ON_IOERROR, CatchIt

IF (Byte(inputValue))[0] EQ 32B THEN RETURN, ""

CASE dataType OF
   'INT': IF inputValue EQ "" OR inputValue EQ "-" OR inputValue EQ "+" THEN $
      retValue = 'NULLVALUE' ELSE retValue = Fix(inputValue)
   'LONG': IF inputValue EQ "" OR inputValue EQ "-" OR inputValue EQ "+" THEN $
      retValue = 'NULLVALUE' ELSE retValue = Long(inputValue)
   'FLOAT' : IF inputValue EQ "" OR inputValue EQ "-" OR inputValue EQ "+" THEN $
      retValue = 'NULLVALUE' ELSE retValue = Float(inputValue)
   'DOUBLE': IF inputValue EQ "" OR inputValue EQ "-" OR inputValue EQ "+" THEN $
      retValue = 'NULLVALUE' ELSE retValue = Double(inputValue)
   'STRING' : retValue = inputValue
ENDCASE
RETURN, retValue

CatchIt:
   retValue = 'NULLVALUE'
   RETURN, retValue
END ;----------------------------------------------------------------------------



Function COYOTE_Field_Validate, value, dataType, Decimal=decimal, Digits=digits, $
   Positive=positive

; This function eliminates illegal characters from a string that represents
; a number. The return value is a properly formatted string that can be turned into
; an INT, LONG, FLOAT, or DOUBLE value.
;
; + 43B
; - 45B
; . 46B
; 0 - 9 48B -57B
; 'eEdD' [101B, 69B, 100B, 68B]

   ; A null string should be returned at once.

IF N_Elements(value) EQ 0 THEN value = ""
value = value[0]
IF value EQ "" THEN RETURN, String(value)

   ; No leading or trainnig blank characters to evaluate.

value = StrTrim(value, 2)

IF N_Elements(dataType) EQ 0 THEN dataType = 'STRING'

   ; A string value should be returned at once. Nothing to check.

IF StrUpCase(dataType) EQ 'STRING' THEN RETURN, value

   ; Check integers and longs. A "-" or "+" in the first character is allowed. Otherwise,
   ; only number between 0 and 9, or 43B to 57B.

IF StrUpCase(dataType) EQ 'INT' OR StrUpCase(dataType) EQ 'LONG' THEN BEGIN

   returnValue = Ptr_New(/Allocate_Heap)
   asBytes = Byte(value)
   IF positive THEN BEGIN
      IF (asBytes[0] EQ 43B) OR $
         (asBytes[0] GE 48B AND asBytes[0] LE 57B) THEN *returnValue = [asBytes[0]]
   ENDIF ELSE BEGIN
      IF (asBytes[0] EQ 45B) OR (asBytes[0] EQ 43B) OR $
         (asBytes[0] GE 48B AND asBytes[0] LE 57B) THEN *returnValue = [asBytes[0]]
   ENDELSE
   length = StrLen(asBytes)
   IF length EQ 1 THEN BEGIN
      IF N_Elements(*returnValue) EQ 0 THEN  *returnValue = [32B] ELSE $
            *returnValue = [asBytes[0]]
   ENDIF ELSE BEGIN
      FOR j=1,length-1 DO BEGIN
         IF (asBytes[j] GE 48B AND asBytes[j] LE 57B) THEN BEGIN
            IF N_Elements(*returnValue) EQ 0 THEN  *returnValue = [asBytes[j]] ELSE $
               *returnValue = [*returnValue, asBytes[j]]
         ENDIF
      ENDFOR
  ENDELSE
  IF N_Elements(*returnValue) NE 0 THEN retValue = String(*returnValue) ELSE retValue = ""
  Ptr_Free, returnValue

      ; Check for digit restrictions.

  IF digits GT 0 THEN BEGIN
      retValue = StrTrim(retValue, 2)
      IF StrMid(retValue, 0, 1) EQ "-" THEN digits = digits + 1
      retValue = StrMid(retValue, 0, digits)
  ENDIF

  RETURN, retValue

ENDIF

   ; Check floating and double values. (+,-) in first character or after 'eEdD'.
   ; Only numbers, signs, decimal points, and 'eEdD' allowed.

IF StrUpCase(dataType) EQ 'FLOAT' OR StrUpCase(dataType) EQ 'DOUBLE' THEN BEGIN
   returnValue = Ptr_New(/Allocate_Heap)
   asBytes = Byte(value)
   IF positive THEN BEGIN
      IF (asBytes[0] EQ 43B) OR $
         (asBytes[0] GE 48B AND asBytes[0] LE 57B) OR $
         (asBytes[0] EQ 46B) THEN *returnValue = [asBytes[0]]
      IF (asBytes[0] EQ 46B) THEN haveDecimal = 1 ELSE haveDecimal = 0
   ENDIF ELSE BEGIN
      IF (asBytes[0] EQ 45B) OR (asBytes[0] EQ 43B) OR $
         (asBytes[0] GE 48B AND asBytes[0] LE 57B) OR $
         (asBytes[0] EQ 46B) THEN *returnValue = [asBytes[0]]
      IF (asBytes[0] EQ 46B) THEN haveDecimal = 1 ELSE haveDecimal = 0
   ENDELSE
   haveExponent = 0
   length = StrLen(asBytes)
   prevByte = asBytes[0]
   exponents = Byte('eEdD')
   IF length EQ 1 THEN BEGIN
      IF N_Elements(*returnValue) EQ 0 THEN  *returnValue = [32B] ELSE $
            *returnValue = [asBytes[0]]
   ENDIF ELSE BEGIN
      FOR j=1,length-1 DO BEGIN
         IF (asBytes[j] GE 48B AND asBytes[j] LE 57B) THEN BEGIN
            IF N_Elements(*returnValue) EQ 0 THEN  *returnValue = [asBytes[j]] ELSE $
               *returnValue = [*returnValue, asBytes[j]]
            prevByte = asBytes[j]
         ENDIF ELSE BEGIN

            ; What kind of thing is it?

            IF (asBytes[j] EQ 46B) THEN BEGIN ; A decimal point.
               IF haveDecimal EQ 0 THEN BEGIN
                  *returnValue = [*returnValue, asBytes[j]]
                  haveDecimal = 1
                  prevByte = asBytes[j]
               ENDIF
            ENDIF

            IF (asBytes[j] EQ 45B) OR (asBytes[j] EQ 43B) THEN BEGIN ; A + or - sign.
               index = Where(exponents EQ prevByte, count)
               IF count EQ 1 AND haveExponent THEN BEGIN
                  *returnValue = [*returnValue, asBytes[j]]
                  aveDecimal = 1
                  prevByte = asBytes[j]
               ENDIF
            ENDIF

            index = Where(exponents EQ asBytes[j], count)
            IF count EQ 1 AND haveExponent EQ 0 THEN BEGIN ; An exponent
               *returnValue = [*returnValue, asBytes[j]]
               haveExponent = 1
               prevByte = asBytes[j]
            ENDIF
         ENDELSE
      ENDFOR
   ENDELSE
   IF N_Elements(*returnValue) NE 0 THEN BEGIN

      retValue = String(*returnValue)
      retValue = StrTrim(retValue, 2)

               ; Check for decimal restrictions

      IF decimal GE 0 THEN BEGIN
         theDecimalPt = StrPos(retValue, '.')
         IF theDecimalPt NE -1 THEN retValue = StrMid(retValue, 0, theDecimalPt + decimal + 1)
      ENDIF

   ENDIF ELSE retValue = ""
   Ptr_Free, returnValue

      ; Is this a representable number?

   testValue = COYOTE_Field_ReturnValue(retValue, dataType)
   IF String(testValue) NE 'NULLVALUE' THEN numCheck = Finite(testValue) ELSE numCheck = 1
   IF numCheck THEN BEGIN
      RETURN, retValue
   ENDIF ELSE BEGIN
      Message, 'The requested number is not representable.', /Informational
      RETURN, ""
   ENDELSE
ENDIF

END ;----------------------------------------------------------------------------



Pro COYOTE_Field__Define

; The COYOTE_Field Event Structure.

   event = { COYOTE_FIELD, $      ; The name of the event structure.
             ID: 0L, $            ; The ID of the compound widget's top-level base.
             TOP: 0L, $           ; The widget ID of the top-level base of the hierarchy.
             HANDLER: 0L, $       ; The event handler ID. Filled out by IDL.
             Value: Ptr_New(), $  ; A pointer to the widget value.
             Type:"" $            ; A string indicating the type of data in the VALUE field.
           }                      ; Values are "INT", "LONG", "FLOAT", "DOUBLE", or "STRING".

END ;----------------------------------------------------------------------------



Pro COYOTE_Field_Kill_Notify, ID

; This routine cleans up the pointer when the compound widget is destroyed.

Widget_Control, ID, Get_UValue=info, /No_Copy
Ptr_Free, info.theValue
END ;----------------------------------------------------------------------------



Pro COYOTE_Field_Set_Value, cw_tlb, value

; This procedure sets a value for the compound widget. The value
; is a value appropriate for the data type or a string.

   ; Get info structure.

info_carrier = Widget_Info(cw_tlb, Find_by_UName='INFO_CARRIER')
Widget_Control, info_carrier, Get_UValue=info, /No_Copy

   ; Validate the value.

theText = Strtrim(value, 2)
theText = COYOTE_Field_Validate(theText, info.dataType, Decimal=info.decimal, $
   Digits=info.digits, Positive=info.positive)

   ; Load the value in the widget.

Widget_Control, info.textID, Set_Value=theText, Set_Text_Select=[StrLen(theText),0]
info.theText = theText

   ; Set the actual value of the compound widget.

*info.theValue = COYOTE_Field_ReturnValue(info.theText, info.dataType)
Widget_Control, info_carrier, Set_UValue=info, /No_Copy
END ;----------------------------------------------------------------------------



Function COYOTE_Field_Get_Value, cw_tlb

; This function returns the numerical or string value of the
; compound widget.

info_carrier = Widget_Info(cw_tlb, Find_by_UName='INFO_CARRIER')
Widget_Control, info_carrier, Get_UValue=info, /No_Copy
value = info.theValue
Widget_Control, info_carrier, Set_UValue=info, /No_Copy
RETURN, value
END ;----------------------------------------------------------------------------




PRO COYOTE_Field_Event_Handler, event

; The main event handler for the compound widget.

   ; Get the info structure. Get the previous text, the current
   ; cursor location in the text widget, and indicate this is not
   ; a Carriage Return event.

Widget_Control, event.ID, Get_UValue=info, /No_Copy
previousText = info.theText
textLocation = Widget_Info(event.id, /Text_Select)
cr_event = 0

   ; What kind of event is this?

possibleTypes = ['INSERT SINGLE CHARACTER', 'INSERT MULTIPLE CHARACTERS', 'DELETE TEXT', 'SELECT TEXT']
thisType = possibleTypes[event.type]

   ; Branch on event type.

CASE thisType OF

   'INSERT SINGLE CHARACTER': BEGIN

            ; Get the current contents of text widget. Validate it.

         Widget_Control, info.textID, Get_Value=newText
         newText = newText[0]
         validText = COYOTE_Field_Validate(newText, info.dataType, Decimal=info.decimal, $
               Digits=info.digits, Positive=info.positive)

            ; If it is valid, leave it alone. If not, go back to previous text.

         IF validText NE newText THEN BEGIN
            Widget_Control, info.textID, Set_Value=previousText, Set_Text_Select=[textLocation[0]-1,0]
         ENDIF ELSE BEGIN
            info.theText = validText
            testValue  = COYOTE_Field_ReturnValue(validText, info.dataType)
            IF String(testValue) EQ "NULLVALUE" THEN BEGIN
               Ptr_Free, info.theValue
               info.theValue = Ptr_New(/Allocate_Heap)
            ENDIF ELSE *info.theValue = testValue
         ENDELSE

            ; Is this a Carriage Return event?

         IF event.ch EQ 10B then cr_event = 1
      ENDCASE

   'INSERT MULTIPLE CHARACTERS': BEGIN

            ; Same thing as above, but for all the characters you are inserting.

         Widget_Control, info.textID, Get_Value=newText
         newText = newText[0]
         validText = COYOTE_Field_Validate(newText, info.dataType, Decimal=info.decimal, $
            Digits=info.digits, Positive=info.positive)
         IF validText NE newText THEN BEGIN
            Widget_Control, info.textID, Set_Value=previousText, Set_Text_Select=[textLocation[0]-1,0]
         ENDIF ELSE BEGIN
            info.theText = validText
            testValue  = COYOTE_Field_ReturnValue(validText, info.dataType)
            IF String(testValue) EQ "NULLVALUE" THEN BEGIN
               Ptr_Free, info.theValue
               info.theValue = Ptr_New(/Allocate_Heap)
            ENDIF ELSE *info.theValue = testValue
         ENDELSE
      ENDCASE

   'DELETE TEXT': BEGIN

            ; Just get the new text and update the info stucture.

         Widget_Control, info.textID, Get_Value=theText
         theText = theText[0]
         validText = COYOTE_Field_Validate(theText, info.dataType, Decimal=info.decimal, $
            Digits=info.digits, Positive=info.positive)

           ; Load the valid text.

        Widget_Control, info.textID, Set_Value=validText, Set_Text_Select=[textLocation[0],0]
        info.theText = validText
        testValue  = COYOTE_Field_ReturnValue(info.theText, info.dataType)
        IF String(testValue) EQ "NULLVALUE" THEN BEGIN
            Ptr_Free, info.theValue
            info.theValue = Ptr_New(/Allocate_Heap)
        ENDIF ELSE *info.theValue = testValue
      ENDCASE

   'SELECT TEXT': ; Nothing to do.

ENDCASE

   ; Do you report all events, or only Carriage Return events?

IF info.cr_only THEN BEGIN
   IF info.event_func NE "" THEN BEGIN
      thisEvent = {COYOTE_Field, info.cw_tlb, event.top, 0L, info.theValue, info.dataType}
      IF cr_event THEN Widget_Control, info.cw_tlb, Send_Event=thisEvent
   ENDIF

   IF info.event_pro NE "" THEN BEGIN
      thisEvent = {COYOTE_Field, info.cw_tlb, event.top, 0L, info.theValue, info.dataType}
      IF cr_event THEN Widget_Control, info.cw_tlb, Send_Event=thisEvent
   ENDIF
ENDIF ELSE BEGIN
   IF info.event_func NE "" THEN BEGIN
      thisEvent = {COYOTE_Field, info.cw_tlb, event.top, 0L, info.theValue, info.dataType}
      Widget_Control, info.cw_tlb, Send_Event=thisEvent
   ENDIF

   IF info.event_pro NE "" THEN BEGIN
      thisEvent = {COYOTE_Field, info.cw_tlb, event.top, 0L, info.theValue, info.dataType}
      Widget_Control, info.cw_tlb, Send_Event=thisEvent
   ENDIF
ENDELSE

   ; Out of here.

Widget_Control, event.ID, Set_UValue=info, /No_Copy
END ;----------------------------------------------------------------------------



Function COYOTE_Field, $          ; The compound widget COYOTE_Field.
   parent, $                      ; The parent widget. Required for all compound widgets.
   Column=column, $               ; Set this keyword to have Label above Text Widget.
   CR_Only=cr_only, $             ; Set this keyword if you only want Carriage Return events.
   Digits=digits, $               ; Set this keyword to number of allowed digits in INT and LONG values.
   Decimal=decimal, $             ; Set to the number of digits to right of decimal point.
   DoubleValue=doublevalue, $     ; Set this keyword if you want DOUBLE values returned.
   Event_Func=event_func, $       ; Set this keyword to the name of an Event Function.
   Event_Pro=event_pro, $         ; Set this keyword to the name of an Event Procedure.
   FieldFont=fieldfont, $         ; The font name for the text in the Text Widget.
   FloatValue=floatvalue, $       ; Set this keyword for FLOAT values.
   Frame=frame, $                 ; Set this keyword to put a frame around the compound widget.
   IntegerValue=integervalue, $   ; Set this keyword for INTEGER values.
   LabelFont=labelfont, $         ; The fon name for the text in the Label Widget.
   LabelSize=labelsize, $         ; The X screen size of the Label Widget.
   LongValue=longvalue, $         ; Set this keyword for LONG values.
   Positive=positive, $           ; Set this keyword to only allow positive number values.
   Row=row, $                     ; Set this keyword to have the Label beside the Text Widget. (The default.)
   Scr_XSize=scr_xsize, $         ; The X screen size of the compound widget.
   Scr_YSize=scr_ysize, $         ; The Y screen size of the compound widget.
   StringValue=stringvalue, $     ; Set this keyword for STRING values. (The default.)
   Title=title, $                 ; The text to go on the Label Widget.
   UValue=uvalue, $               ; A user value for any purpose.
   Value=value, $                 ; The "value" of the compound widget.
   XSize=xsize                    ; The X size of the Text Widget.

   ; A parent is required.

IF N_Elements(parent) EQ 0 THEN BEGIN
   Message, 'A PARENT argument is required. Returning...', /Informational
   RETURN, -1L
ENDIF

   ; Check keyword values.

IF N_Elements(column) EQ 0 THEN column = 0
IF N_Elements(digits) EQ 0 THEN digits = 0 ELSE digits = Fix(digits)
IF N_Elements(decimal) EQ 0 THEN decimal = -1 ELSE decimal = Fix(decimal)
IF N_Elements(event_func) EQ 0 THEN event_func = ""
IF N_Elements(event_pro) EQ 0 THEN event_pro = ""
IF N_Elements(fieldfont) EQ 0 THEN fieldfont = ""
IF N_Elements(frame) EQ 0 THEN frame = 0
IF N_Elements(labelfont) EQ 0 THEN labelfont = ""
IF N_Elements(labelsize) EQ 0 THEN labelsize = 0
IF N_Elements(scr_xsize) EQ 0 THEN scr_xsize = 0
IF N_Elements(scr_ysize) EQ 0 THEN scr_ysize = 0
IF N_Elements(title) EQ 0 THEN title = "Input Value: "
IF N_Elements(uvalue) EQ 0 THEN uvalue = ""
IF N_Elements(value) EQ 0 THEN value = "" ELSE value = StrTrim(value,2)
IF N_Elements(xsize) EQ 0 THEN xsize = 0
IF N_Elements(row) EQ 0 AND column EQ 0 THEN row = 1 ELSE row = 0
positive = Keyword_Set(positive)

   ; What data type are we looking for?

dataType = 'STRING'
IF Keyword_Set(stringvalue) THEN dataType = 'STRING'
IF Keyword_Set(integervalue) THEN dataType = 'INT'
IF Keyword_Set(longvalue) THEN dataType = 'LONG'
IF Keyword_Set(floatvalue) THEN dataType = 'FLOAT'
IF Keyword_Set(doublevalue) THEN dataType = 'DOUBLE'

   ; Validate the input value.

value = COYOTE_Field_Validate(value, dataType, Decimal=decimal, Digits=digits, Positive=positive)

   ; Create the widgets.

cw_tlb = Widget_Base( parent, $               ; The top-level base of the compound widget.
   Pro_Set_Value='COYOTE_Field_Set_Value', $
   Func_Get_Value='COYOTE_Field_Get_Value', $
   Frame=frame, $
   Row=row, $
   Column=Keyword_Set(column), $
   Base_Align_Center=1, $
   UValue=uvalue, $
   Event_Pro=event_pro, $
   Event_Func=event_func )

labelID = Widget_Label( cw_tlb, Value=title, Font=labelfont, $ ; The Label Widget.
  Scr_XSize=labelsize)

textID = Widget_Text( cw_tlb, $  ; The Text Widget.
   Value=value, $
   XSize=xsize, $
   YSize=1, $
   Scr_XSize=scr_xsize, $
   Scr_YSize=scr_ysize, $
   Font=fieldfont, $
   All_Events=1, $
   Event_Pro='COYOTE_Field_Event_Handler', $
   UName='INFO_CARRIER', $
   Kill_Notify='COYOTE_Field_Kill_Notify', $
   Editable=1 )

   ; Set the actual return value of the compound widget.

theValue = COYOTE_Field_ReturnValue(value, dataType)

   ; The info structure.

info = { theText:value, $                  ; The text in the Text Widget.
         theValue:Ptr_New(theValue), $     ; The real value of the Text Widget.
         cw_tlb:cw_tlb, $                  ; The top-level base of the compound widget.
         event_func:event_func, $          ; The name of the event handler function.
         event_pro:event_pro, $            ; The name of the event handler procedure.
         cr_only:Keyword_Set(cr_only), $   ; A flag to return events only on CR events.
         dataType:dataType, $              ; The type of data wanted in the compound widget.
         decimal:decimal, $                ; The number of digits in decimal numbers.
         digits:digits, $                  ; The number of digits in integer numbers.
         positive:positive, $              ; Flag to indicate positive number values.
         textID:textID }                   ; The widget identifier of the Text Widget.

    ; Store info structure in Text Widget.

Widget_Control, textID, Set_UValue=info, /No_Copy
RETURN, cw_tlb
END ;----------------------------------------------------------------------------



PRO GETIMAGE_CenterTLB, tlb

Device, Get_Screen_Size=screenSize
xCenter = screenSize(0) / 2
yCenter = screenSize(1) / 2

geom = Widget_Info(tlb, /Geometry)
xHalfSize = geom.Scr_XSize / 2
yHalfSize = geom.Scr_YSize / 2

Widget_Control, tlb, XOffset = xCenter-xHalfSize, $
   YOffset = yCenter-yHalfSize

END ;-------------------------------------------------------------------------



PRO GETIMAGE_NULL_EVENTS, event

   ; The purpose of this event handler is to do nothing
   ;and ignore all events that come to it.

END ;-------------------------------------------------------------------------



FUNCTION GETIMAGE_FIND_COYOTE

   ; The purpose of this function is to find the "coyote"
   ; training directory and return its path. If no
   ; directory is found, the function returns a null string.

ON_ERROR, 1

   ; Check this directory first.

CD, Current=thisDir
IF STRPOS(STRUPCASE(thisDir), 'COYOTE') GT 0 THEN RETURN, thisDir

   ; Look in !Path directories.

pathDir = EXPAND_PATH(!Path, /Array)
s = SIZE(pathDir)
IF s(1) LT 1 THEN RETURN, ''
FOR j=0,s(1)-1 DO BEGIN
   check = STRPOS(STRUPCASE(pathDir(j)), 'COYOTE')
   IF check GT 0 THEN RETURN, pathDir(j)
ENDFOR
RETURN, ''
END ;-------------------------------------------------------------------------



PRO GETIMAGE_EVENT, event

   ; The only events that can come here are button events.

   ; Get the info structure out of the user value of the top-level base.

Widget_Control, event.top, Get_UValue=info

   ; There may be errors we can't anticipate. Catch them here, alert the
   ; user as to what the error was, and exit the event handler without
   ; doing any damage.

Catch, error
IF error NE 0 THEN BEGIN
   ok = Widget_Message(!Err_String)
   formdata = {cancel:1}
   *info.ptrToFormData =formdata
   Widget_Control, event.top, /Destroy
   RETURN
ENDIF

   ; Which button caused this event?

Widget_Control, event.id, Get_Value=buttonValue

CASE buttonValue OF

   'Pick Filename': BEGIN

         ; Start in the directory listed in the directory text widget.
         ; Convert the text value to a scalar.

      Widget_Control, info.dirnameID, Get_Value=startDirectory
      startDirectory = startDirectory(0)

         ; If this directory doesn't exist, use the current directory.

      test = Findfile(startDirectory, Count=foundfile)
      IF foundfile NE 1 THEN CD, Current=startDirectory

         ; Use PICKFILE to pick a name.

      pick = Dialog_Pickfile(Path=startDirectory, /NoConfirm, $
         Get_Path=path, Filter='*.')

         ; Set the directory text widget with the name of the directory.
         ; Make sure the user didn't cancel out of PICKFILE.

      IF pick NE '' THEN BEGIN

            ; Find the lengths of the PICK and the PATH.

         pathLen = StrLen(path)
         picklen = StrLen(pick)

           ; Shorten the PATH to take off last file separator.

         path = StrMid(path,0,pathLen-1)

            ; Put the PATH in the directory location.

         Widget_Control, info.dirnameID, Set_Value=path

            ; Set the filename text widget with the name of the file.

         filename = StrMid(pick, pathlen, picklen-pathlen)
         Widget_Control, info.filenameID, Set_Value=filename

      ENDIF

      END ; of the Pick Filename button case

    'Cancel': BEGIN

         ; Have to exit here gracefully. Set the "CANCEL" flag.

      formdata = {cancel:1}
      *info.ptrToFormData =formdata

         ; Out of here!

      Widget_Control, event.top, /Destroy
      END ; of the Cancel button case

    'Accept': BEGIN  ; Gather the form information.

          ; Put the directory and filename together to make a path.

       Widget_Control, info.dirnameID, Get_Value=directory
       Widget_Control, info.filenameID, Get_Value=file

       filename = Filepath(Root_Dir=directory(0),file(0))

          ; Get the size and header info. Remember these are Pointers to Integers!

       Widget_Control, info.headerID, Get_Value=header
       Widget_Control, info.xsizeID, Get_Value=xsize
       Widget_Control, info.ysizeID, Get_Value=ysize
       Widget_Control, info.frameID, Get_Value=frames

       header = Fix(*header)
       xsize =  Fix(*xsize)
       ysize =  Fix(*ysize)
       frames =  Fix(*frames)

          ; Get the data type from the droplist widget.

       listIndex = Widget_Info(info.droplistID, /Droplist_Select)
       datatype = info.datatypes(listIndex)

          ; Get the format index from the formatlist widget.

       formatIndex = Widget_Info(info.formatlistID, /Droplist_Select)

          ; Create the formdata structure from the information you collected.

       formdata = {header:header, xsize:xsize, ysize:ysize, frames:frames, $
          filename:filename, datatype:datatype, formatIndex:formatIndex, cancel:0}

          ; Store the formdata in the pointer location.

       *info.ptrToFormData = formdata

         ; Out of here!

      Widget_Control, event.top, /Destroy
      END ; of the Accept button case

ENDCASE

END ; of GETIMAGE_EVENT event handler ***************************************



FUNCTION GETIMAGE, filename, Directory=directory, XSize=xsize, YSize=ysize, $
   Frames=frames, Header=header, Parent=parent, XDR=xdr, XOffSet=xoffset, $
   YOffSet=yoffset, Cancel=canceled, Catch=catchit, Group_Leader=group_leader

   ; This is a function to specify the size, data type, and header information
   ; about an image that you would like to read. It reads the data and returns
   ; it as the result of the function. If an error occurs or the user CANCELS,
   ; the function returns a 0.

   ; Must have IDL 5 because of pointers and other functionality.

On_Error, 2

thisRelease = StrMid(!Version.Release, 0, 1)
IF thisRelease NE '5' THEN BEGIN
   ok = Widget_Message('This program requires IDL 5 functionality. Sorry.')
   RETURN, 0
ENDIF

   ; Catch errors unless explicitly told not to.

IF N_Elements(catchit) EQ 0 THEN catchit = 1

   ; Respond to either PARENT or GROUP_LEADER keywords.

IF N_Elements(parent) NE 0 THEN group_leader = parent

   ; Check for parameters and keywords.

IF N_Params() EQ 0 THEN filename='ctscan.dat'

   ; If DIRECTORY keyword is not used, use the "coyote" directory.
   ; If that is not found, use the current directory.

IF N_Elements(directory) EQ 0 THEN BEGIN

   startDirectory = GetImage_Find_Coyote()
   IF startDirectory EQ '' THEN BEGIN
      ;CD, Current=startDirectory
      dir = Filepath(Subdirectory=['examples', 'data'], 'ctscan.dat')
      startDirectory = StrMid(dir, 0, StrLen(dir)-11)
   ENDIF

ENDIF ELSE startDirectory = directory

   ; If the default file is not in the directory, make the filename
   ; a null string.

thisFile = Filepath(Root_Dir=startDirectory, filename)
ok = Findfile(thisFile, Count=count)
IF count EQ 0 THEN filename = ''

   ; Check for size and header keywords. These probably come in as
   ; numbers and you need strings to put them into text widgets.

IF N_Elements(xsize) EQ 0 THEN xsize='256' ELSE xsize=StrTrim(xsize,2)
IF N_Elements(ysize) EQ 0 THEN ysize='256' ELSE ysize=StrTrim(ysize,2)
IF N_Elements(frames) EQ 0 THEN frames='0' ELSE frames=StrTrim(frames,2)
IF N_Elements(header) EQ 0 THEN header='0' ELSE header=StrTrim(header,2)

   ; Create a modal top-level base if group leader is present.

IF N_Elements(group_leader) EQ 0 THEN $
   tlb = Widget_Base(Column=1, Title='Read Image Data', /Base_Align_Center) ELSE $
   tlb = Widget_Base(Column=1, Title='Read Image Data', Modal=1, $
      Group_Leader=group_leader, /Base_Align_Center)

frameBase = Widget_Base(tlb, Frame=1, Column=1)

   ; Create the directory widgets.

dirnamebase = Widget_Base(frameBase, Row=1)
   dirnamelabel = Widget_Label(dirnamebase, Value='Directory:')
   dirnameID = Widget_Text(dirnamebase, Value=startDirectory, /Editable, $
      Event_Pro='GETIMAGE_NULL_EVENTS', XSize=Fix(2.0*StrLen(startDirectory) > 50))

   ; Create the filename widgets.

filenamebase = Widget_Base(frameBase, Row=1)
   filenamelabel = Widget_Label(filenamebase, Value='Filename:')
   filenameID = Widget_Text(filenamebase, Value=filename, /Editable, $
      Event_Pro='GETIMAGE_NULL_EVENTS', XSize=2*StrLen(filename) > 20)

   ; Create a button to allow user to pick a filename.

pickbutton = Widget_Button(filenamebase, Value='Pick Filename')

   ; Create a droplist widget to select file data types.

database = Widget_Base(frameBase, Row=1)
   datatypes = ['Byte', 'Integer', 'Long', 'Float']
   droplistID = Widget_Droplist(database, Value=datatypes, $
      Title='Data Type: ', Event_Pro='GETIMAGE_NULL_EVENTS')

   ; Create a droplist widget to select file formats.

   formatlistID = Widget_Droplist(database, Value=['None', 'XDR'], $
      Title='File Format: ', Event_Pro='GETIMAGE_NULL_EVENTS')

   ; Create a text widget to accept a header size.

   headerID = Coyote_Field(database, Value=header, XSize=8, Title='Header Size:', /Integer)

   ; Create widgets to gather the required file sizes.

sizebase = Widget_Base(frameBase, Row=1)

   xsizeID = Coyote_Field(sizebase, Title='X Size:', Value=xsize, /Integer, XSize=8, /Positive)

   ysizeID = Coyote_Field(sizebase, Title='Y Size:', Value=ysize, /Integer, XSize=8, /Positive)

   frameID = Coyote_Field(sizebase, Title='Frames:', Value=frames, /Integer, XSize=8, /Positive)

   ; Create cancel and accept buttons.

cancelbase = Widget_Base(tlb, Row=1)
   cancel = Widget_Button(cancelbase, Value='Cancel')
   accept = Widget_Button(cancelbase, Value='Accept')

   ; Recalculate the length of the filenameID widget.

g1 = Widget_Info(dirnamebase, /Geometry)
g2 = Widget_Info(filenamelabel, /Geometry)
g3 = Widget_Info(pickbutton, /Geometry)
target = g1.scr_xsize - (g2.scr_xsize + g3.scr_xsize) - 10
Widget_Control, filenameID, Scr_XSize=target

   ; Center the program on the display.

GetImage_CenterTLB, tlb

   ; Realize the program.

Widget_Control, tlb, /Realize

   ; Create a pointer to store the information collected from the form.
   ; The initial data stored here is set to CANCEL, so nothing needs to
   ; be done if the user kills the widget with the mouse..

ptrToFormData = Ptr_New({cancel:1})

   ; Set the correct file format in the format droplist widget.

Widget_Control, formatlistID, Set_Droplist_Select=Keyword_Set(xdr)

   ; Set the text insertion point at the end of the filename text widget.

tip = [StrLen(filename),0]
Widget_Control, filenameID, Input_Focus=1
Widget_Control, filenameID, Set_Text_Select=tip

   ; Create an info structure with program information.

info = { filenameID:filenameID, $        ; The name of the file.
         dirnameID:dirnameID, $          ; The ID of the widget with the directory name.
         xsizeID:xsizeID, $              ; The ID of the widget with the image X size.
         ysizeID:ysizeID, $              ; The ID of the widget with the image Y size.
         frameID:frameID, $              ; The ID of the widget with the image FRAME size.
         headerID:headerID, $            ; The ID of the widget with the image HEADER size.
         droplistID:droplistID, $        ; The ID of the image data TYPE droplist ID.
         formatlistID:formatlistID, $    ; The ID of the image FORMAT droplist ID.
         datatypes:datatypes, $          ; The possible data types.
         ptrToFormData:ptrToFormData}    ; A pointer to store the form information.

  ; Store the info structure in the user value of the top-level base.

Widget_Control, tlb, Set_UValue=info

   ; The form will be a MODAL or BLOCKING widget, depending upon the
   ; presence of the PARENT. We block or are MODAL here until the widget
   ; is destroyed.

XManager, 'getimage', tlb, Event_Handler='GETIMAGE_EVENT'

   ; When the widget is destroyed, the block is released, and we
   ; return here. Get the form data that was collected by the form
   ; and stored in the pointer location.

formdata = *ptrToFormData

   ; If there is nothing here. Free the pointer and return.

IF N_Elements(formdata) EQ 0 THEN BEGIN
   Ptr_Free, ptrToFormData
   canceled = 1
   RETURN, 0
ENDIF

   ; Did the user cancel out of the form? If so, return a 0.

IF formdata.cancel EQ 1 THEN BEGIN
   Ptr_Free, ptrToFormData
   canceled = 1
   RETURN, 0
ENDIF

   ; Make the proper sized image array. Check for success.

image = 0
IF STRUPCASE(formdata.datatype) EQ 'INTEGER' THEN formdata.datatype = 'INT'
IF formdata.frames EQ 0 THEN $
   command = 'image = Make_Array(formdata.xsize, formdata.ysize, ' + $
      formdata.datatype + '=1)' ELSE $
   command = 'image = Make_Array(formdata.xsize, ' + $
      'formdata.ysize, formdata.frames, ' + formdata.datatype + '=1)'

check = Execute(command)
IF check EQ 0 THEN BEGIN
   ok = Dialog_Message("Problem making image array. Returning 0.")
   canceled = 1
   RETURN, 0
ENDIF

   ; We can have all kinds of trouble reading data. Let's catch all
   ; input and output errors and alert user without crashing the program!

IF Keyword_Set(catchit) THEN Catch, theError ELSE theError = 0
IF theError NE 0 THEN BEGIN

      ; If we can't read the file for some reason, let the user know
      ; why, free the pointer and its information, check the logical
      ; unit number back in if it is checked out, and return a 0.

   ok = Dialog_Message(!Err_String)
   Ptr_Free, ptrToFormData
   canceled = 1
   IF N_ELements(lun) NE 0 THEN Free_Lun, lun
   RETURN, 0
ENDIF

   ; Set the canceled flag.

canceled = formdata.cancel

   ; Read the data file.

IF formdata.header GT 0 THEN header = BytArr(formdata.header)
Get_Lun, lun
OpenR, lun, formdata.filename, XDR=formdata.formatIndex
IF formdata.header EQ 0 THEN ReadU, lun, image $
   ELSE ReadU, lun, header, image
Free_Lun, lun

   ; Free the pointer.

Ptr_Free, ptrToFormData

RETURN, image

END ; of GETIMAGE program ***************************************************
