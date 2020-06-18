;+
; NAME:
;	CW_ImageFull
; PURPOSE:
;	Defines a compound widget containing a draw widget to be
;	used to display entire images at reduced resolution.
; CALLING SEQUENCE:
;	cw = CW_ImageFull(Img [, Parent])
; INPUTS:
;	Img    - The image to be displayed.
;	Parent - The parent widget.  If not specified, this is a top level
;	         widget.
; RETURNED:
;	cw     - The widget identifier for the compound widget.
; KEYWORDS:
;	Button_Events - If present and nonzero, mouse button events are
;	                passed back to the application and/or parent
;	                widget.
;	BytScl        - If present and nonzero, the image is bytescaled
;	                before it is displayed.
;	Retain        - Allows the user to specify the backing store to
;	                use.  Defaults to 2 -- IDL provides backing store.
;	ShowBox       - If present and nonzero, support for drawing a box
;	                in the draw widget is supplied.
;	X_Size        - Specifies the horizontal size of the draw widget
;	                display in pixels.  Defaults to 512.
;	Y_Size        - Specifies the vertical size of the draw widget
;	                display in pixels.  Defaults to 512.
;	_EXTRA        - Any allowed Widget_Base keywords may be supplied
;	                through IDL keyword inheritance; these keywords
;	                are passed along to the Widget_Base that provides
;	                the basis of the compound widget.
; COMMENTS:
;	This very simple compound widget simply wraps around a draw widget
;	displaying an image to be displayed.  The image is resized to fit
;	within the draw widget using CONGRID.
;
;	If the BytScl keyword is specified, the image is loaded into the
;	draw widget graphics window using TVSCL; else TV is used.
;
;	Draw widget events are passed back in a structure CW_IMAGEFULL,
;	containing:
;		Id      - The id of the compound widget.
;		Top     - The id of the 
;		Handler - The id of the widget owning the event handler.
;		          Always set to 0 here, but changed by the IDL
;		          event handling system.
;		Type    - The type of event: 0=mouse button press, 1=mouse
;		          button release, 2=mouse motion, 3=scrollbars moved,
;		          4=visibility changed.  Note: only 0 and 1 are
;		          currently enabled to be passed back.
;		X       - The horizontal image position of the event.
;		Y       - The vertical image position of the event.
;		Press   - Identifies the button(s) pressed, with the
;		          left button represented by the least sig. bit.
;		Release - Identifies the button(s) released, with the
;		          left button represented by the least sig. bit.
;		Clicks  - Double click? 1=no, 2=yes.
;		Value   - The value of the image at the mouse button event.
;		          Returned in double precision regardless of the
;		          data type of the image.
;
;	A box can be drawn on the draw widget using the CW_IMAGEFULL_BOX
;	procedure.  The syntax is:
;	  CW_ImageFull_Box, id, xLL, yLL, xUR, yUR [, /Device]
;	where id is the compound widget id and xLL ... define the corners of
;	the box as indices into the FULL RESOLUTION IMAGE ARRAY.  If Device
;	is specified, the coordinates are assumed to be in device coordinates.
;	The box actually consists of two boxes: the outer one in black (color
;	index 0) and the inner one in white (color index 255).
;
;	As with most compound widgets, it is not recommended that this
;	widget be used as a top level widget.
; MODIFICATION HISTORY:
;	Michael R. Greason, Raytheon ITSS, 03 March 1999
;-
; =============================================================================
;
Pro CW_ImageFull_Box, id, xLL, yLL, xUR, yUR, Device=dev
;
on_error, 2
;
;			Check arguments.
;
If (n_params() LT 5) Then message, $
	'CW_ImageFull_Box, id, xLL, yLL, xUR, yUR [, /Device]'
;
;                       Get the state data.
;
stash = Widget_Info(id, /Child)
Widget_Control, stash, Get_UValue=state, /No_Copy
If ((state.PixMap GE 0) AND state.ShowBox) Then Begin
;
;			make the draw widget the draw widget the current window.
;
	Widget_Control, state.DrawId, Get_Value=win
	oldwin = !d.window
	wset, win
;
;			Determine the corners of the box.
;
;				Device coordinates.
;
	x = [xLL, xLL, xUR, xUR, xLL]
	y = [yLL, yUR, yUR, yLL, yLL]
;
;				Image indices.
;
	If (NOT keyword_set(dev)) Then Begin
		If (state.X_Scale NE 0.0) Then x = x / state.X_Scale
		If (state.Y_Scale NE 0.0) Then y = y / state.Y_Scale
	EndIf
;
;			Reload the base image into the draw widget.
;
	device, Copy=[0, 0, state.X_Win_Size, state.Y_Win_Size, 0, 0, state.PixMap]
;
;			Draw the box.
;
	dx = [+1, +1, -1, -1, +1]
	dy = [+1, -1, -1, +1, +1]
	plots, x, y, Color=0, /Device
	plots, (x + dx), (y + dy), Color=255, /Device
;
;			Reset the current graphics device to the old one.
;
	If (oldwin GE 0) Then wset, oldwin
;
;                       Reset the state data.
;
EndIf
Widget_Control, stash, Set_UValue=state, /No_Copy
;
Return
End
;
; =============================================================================
; Assigns the image displayed in the compound widget.  The local copy of the
; image is resized to the display size of the draw widget using congrid.
; This routine does NOT fill the draw widget.
;
Pro CW_ImageFull_Assign, state, val
;
sz = size(val)
If (sz[0] NE 2) Then Begin
	img = 0B 
	state.X_Scale = 1.0
	state.Y_Scale = 1.0
EndIf Else Begin
	img = congrid(val, state.X_Win_Size, state.Y_Win_Size)
	state.X_Scale = float(sz[1]) / float(state.X_Win_Size)
	state.Y_Scale = float(sz[2]) / float(state.Y_Win_Size)
EndElse
Widget_Control, state.DrawId, Set_UValue=img, /No_Copy
;
Return
End
;
; =============================================================================
; Paints the draw widget.
;
Pro CW_ImageFull_DrawImage, state
;
;			Make the draw widget the current graphics device.
;
Widget_Control, state.DrawId, Get_Value=win
oldwin = !d.window
wset, win
;
;			If the image is 2D, fill the widget with it,
;			resizing the draw widget.
;
Widget_Control, state.DrawId, Get_Uvalue=img, /No_Copy
;
sz = size(img)
If (sz[0] EQ 2) Then Begin
	Widget_Control, state.DrawId, Draw_XSize=sz[1], Draw_YSize=sz[2]
	If (state.Scale) Then Begin
		tvscl, img
	EndIf Else Begin
		tv,    img
	EndElse
	If (state.PixMap GE 0) Then wdelete, state.PixMap
	If (state.ShowBox) Then Begin
		window, /Free, /PixMap, Retain=2, $
			XSize=state.X_Win_Size, YSize=state.Y_Win_Size
		device, Copy=[0, 0, state.X_Win_Size, state.Y_Win_Size, 0, 0, win]
		state.PixMap = !d.Window
		wset, win
	EndIf
EndIf Else Begin
;
;			Otherwise empty the widget.
;
	erase
	If ((state.PixMap GE 0) AND state.ShowBox) Then erase, Channel=state.PixMap
EndElse
;
Widget_Control, state.DrawId, Set_Uvalue=img, /No_Copy
;
;			Reset the current graphics device to the old one.
;
If (oldwin GE 0) Then wset, oldwin
;
Return
End
;
; =============================================================================
; Returns the image value at a given position.
;
Function CW_ImageFull_Value, state, x, y
;
Widget_Control, state.DrawId, Get_Uvalue=img, /No_Copy
;
sz = size(img)
If (sz[0] EQ 2) Then Begin
	xp = (x > 0) < (sz[1] - 1)
	yp = (y > 0) < (sz[2] - 1)
	val = img[xp, yp]
EndIf Else Begin
	val = 0
EndElse
;
Widget_Control, state.DrawId, Set_Uvalue=img, /No_Copy
;
Return, val
End
;
; =============================================================================
; The Widget Realization procedure.  This is called once when the widget
; is realized.  At this time the draw widget has a graphics window assigned
; to it and can be (and is) filled.
;
Pro CW_ImageFull_Notify, id
;
;                       Get the state data.
;
stash = Widget_Info(id, /Child)
Widget_Control, stash, Get_UValue=state, /No_Copy
;
;                       Paint the draw widget.
;
CW_ImageFull_DrawImage, state
;
;                       Reset the state data.
;
Widget_Control, stash, Set_UValue=state, /No_Copy
;
Return
End
;
; =============================================================================
; Returns the image displayed in the compound widget.
;
Function CW_ImageFull_GetValue, id
;
;                       Get the state data.
;
stash = Widget_Info(id, /Child)
Widget_Control, stash, Get_UValue=state, /No_Copy
;
;                       Get the value.
;
Widget_Control, state.DrawId, Get_Uvalue=val
;
;                       Reset the state data.
;
Widget_Control, stash, Set_UValue=state, /No_Copy
;
Return, val
End
;
; =============================================================================
; Assigns the image displayed in the compound widget.  The local copy of the
; image is resized to the display size of the draw widget using congrid.
;
Pro CW_ImageFull_SetValue, id, val
;
;                       Get the state data.
;
stash = Widget_Info(id, /Child)
Widget_Control, stash, Get_UValue=state, /No_Copy
;
;                       Set the value.
;
CW_ImageFull_Assign, state, val
CW_ImageFull_DrawImage, state
;
;                       Reset the state data.
;
Widget_Control, stash, Set_UValue=state, /No_Copy
;
Return
End
;
; =============================================================================
;
Function CW_ImageFull_Event, event
;
ret = 0
tn  = strupcase(strtrim(Tag_Names(event, /Structure_Name), 2))
;
;                      Get the state data.
;
id = event.handler
stash = Widget_Info(id, /Child)
Widget_Control, stash, Get_UValue=state, /No_Copy
;
;		Deal with the event.
;
resflg = 1
Case (tn) Of
;
;			Draw widget events.
;
	'WIDGET_DRAW' : Begin
;
;				Return a mouse position event to the
;				parent or application.
;
		typ = ((event.type EQ 0) OR (event.type EQ 1))
		If (state.BEvent AND typ) Then Begin
		    ret = {CW_IMAGEFULL,	$
			Id:      id,		$
			Top:     event.Top,	$
			Handler: 0L,		$
			Type:    event.Type,	$
			X:       event.X,	$
			Y:       event.Y,	$
			Press:   event.Press,	$
			Release: event.Release,	$
			Clicks:  event.Clicks,	$
			Value:   double(CW_ImageFull_Value(state, event.X, event.Y))}
		    If (state.X_Scale NE 0.0) Then $
		   	ret.X = ret.X * state.X_Scale
		    If (state.Y_Scale NE 0.0) Then $
		   	ret.Y = ret.Y * state.Y_Scale
		EndIf
	End
;
	Else : ret = 0
EndCase
;
;                      Reset the state data.
;
If (resflg NE 0) Then Widget_Control, stash, Set_UValue=state, /No_Copy
;
Return, ret
End
;
; =============================================================================
; Defines the compound widget.
;
Function CW_ImageFull, Img, Parent, BytScl=bscl, 	$
	Button_Events=bevt, Retain=ret, 		$
	X_Size=xscrsz, Y_Size=yscrsz,			$
	ShowBox=sb, _EXTRA=ext
;
on_error, 2
;
;			Check arguments.
;
If (n_elements(xscrsz) LE 0) Then xscrsz = 512
If (n_elements(yscrsz) LE 0) Then yscrsz = 512
If (n_elements(ret)    LE 0) Then ret    =   2
;
sz = size(Img)
If (sz[0] EQ 2) Then Begin
	val = Img
EndIf Else Begin
	val = 0B
EndElse
;
;			Create the widget.
;
evfn  = 'CW_ImageFull_Event'
guvfn = 'CW_ImageFull_GetValue'
suvpr = 'CW_ImageFull_SetValue'
nrpr  = 'CW_ImageFull_Notify'
;
If (n_params() LT 2) Then Begin
	base = Widget_Base(/Column,		$
			Event_Func=evfn,	$
			Func_Get_Value=guvfn,	$
			Pro_Set_Value=suvpr,	$
			Notify_Realize=nrpr,	$
			_EXTRA=ext)
EndIf Else Begin
	base = Widget_Base(Parent, /Column,	$
			Event_Func=evfn,	$
			Func_Get_Value=guvfn,	$
			Pro_Set_Value=suvpr,	$
			Notify_Realize=nrpr,	$
			_EXTRA=ext)
EndElse
;
;				Dummy widgets allowing the draw widget
;				user value to be used.
;
dum = Widget_Base(base)
;
;				The draw widget.
;
drw = Widget_Draw(base,	$
	Retain=ret,	$
	XSize=xscrsz,	$
	YSize=yscrsz,	$
	Button_Events=keyword_set(bevt))
;
;			Define the compound widget internal state.
;
state = {DrawId:  drw,			$
         Scale:   keyword_set(bscl),	$
         BEvent:  keyword_set(bevt),	$
	 X_Win_Size: xscrsz,		$
	 Y_Win_Size: yscrsz,		$
	 X_Scale: 1.0,			$
	 Y_Scale: 1.0,			$
	 ShowBox: keyword_set(sb),	$
	 PixMap:  -1}
CW_ImageFull_Assign, state, val
Widget_Control, Widget_Info(base, /Child), Set_UValue=state, /No_Copy
;
Return, base
End
