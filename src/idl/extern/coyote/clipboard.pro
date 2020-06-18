;+
; NAME:
;       CLIPBOARD
;
; PURPOSE:
;
;       The purpose of this program is to copy the contents of a
;       graphics window to the clipboard for subsequent pasting into
;       applications such as Photoshop or Powerpoint.
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
;      Graphics.
;
; CALLING SEQUENCE:
;
;      CLIPBOARD, window_index
;
; OPTIONAL INPUTS:
;
;       window_index:    The window index number of the graphics window to
;                        copy. If absent, the current graphics window is used
;                        by default.
;
; OUTPUTS:
;       None.
;
; COMMON BLOCKS:
;       None.
;
; DEPENDENCIES:
;
;       Uses the IDLgrClipboard object introduced in IDL 5.2(?).
;
; PROCEDURE:
;
;       Copies the window contents to a clipboard object.
;
; EXAMPLE:
;
;        IDL> Window
;        IDL> Plot, Findgen(11)
;        IDL> CLIPBOARD
;
; RESTRICTIONS:
;
;       May not work for all applications. Applications tested successfully
;       include: Framemaker, Powerpoint, Photoshop, Excel, Microsoft Word.
;       Converts 24-bit images to 2D images with color tables.
;
; MODIFICATION HISTORY:
;
;       Written by: David W. Fanning, 24 October 2001.
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


PRO Clipboard, windowIndex

; This procedure copies the window identified by the
; window index number (or the current window if an index
; number is not provided) to the clipboard.

IF N_Elements(windowIndex) EQ 0 THEN windowIndex = !D.Window

   ; Is this a valid window?

IF windowIndex LT 0 THEN BEGIN
   ok = Dialog_Message('No current window to copy. Returning...')
   RETURN
ENDIF

   ; Catch window setting errors.

Catch, error
IF error NE 0 THEN BEGIN
   Catch, /Cancel
   ok = Dialog_Message('Specified window is unavailable: ' + $
      StrTrim(windowIndex, 2) + '. Returning...')
   WSet, thisWindow
   RETURN
ENDIF

   ; Set active window.

thisWindow = !D.Window
WSet, windowIndex
Catch, /Cancel

   ; Take a snapshot of window. Pay attention to visual depth.

Device, Get_Visual_Depth=thisDepth
IF thisDepth GT 8 THEN BEGIN
   snapshot = TVRD(True=1)
   snapshot = Color_Quan(snapshot, 1, r, g, b)
ENDIF ELSE BEGIN
   snapshot = TVRD()
   TVLCT, r, g, b, /Get
ENDELSE
s = Size(snapshot, /Dimensions)

   ; Create an object graphics image and hierarchy.

palette = Obj_New('IDLgrPalette', r, g, b)
image = Obj_New('IDLgrImage', snapshot, Palette=palette)
model = Obj_New('IDLgrModel')
model->Add, image
thisView = Obj_New('IDLgrView', ViewPlane_Rect=[0,0,s[0],s[1]])
thisView->Add, model

   ; Create a clipboard

theClipboard = Obj_New('IDLgrClipboard', Color_Model=1, $
   Dimensions=[s[0], s[1]], N_Colors=!D.Table_Size, $
   Resolution=[1.0/!D.X_PX_CM, 1.0/!D.Y_PX_CM], $
   Palette=palette)

   ; Copy the snapshot to the clipboard.

theClipboard->Draw, thisView

   ; Destroy the objects.

Obj_Destroy, palette
Obj_Destroy, model
Obj_Destroy, thisView
Obj_Destroy, theClipboard

   ; Restore the current window.

IF thisWindow NE -1 THEN WSet, thisWindow
END