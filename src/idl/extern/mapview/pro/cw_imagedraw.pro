;+
; NAME:
;	CW_ImageDraw
; PURPOSE:
;	Defines a compound widget containing a scrollable draw widget to be
;	used to display images.
; CALLING SEQUENCE:
;	cw = CW_ImageDraw(Img [, Parent])
; INPUTS:
;	Img    - The image to be displayed.
;	Parent - The parent widget.  If not specified, this is a top level
;	         widget.
; RETURNED:
;	cw     - The widget identifier for the compound widget.
; KEYWORDS:
;	Button_Center - If present and nonzero, the view is centered on
;	                to the position of a mouse button click.
;	Button_Events - If present and nonzero, mouse button events are
;	                passed back to the application and/or parent
;	                widget.
;	BytScl        - If present and nonzero, the image is bytescaled
;	                before it is displayed.
;	Retain        - Allows the user to specify the backing store to
;	                use.  Defaults to 2 -- IDL provides backing store.
;	X_Scroll_Size - Specifies the horizontal size of the draw widget
;	                display in pixels.  Defaults to 1024.
;	Y_Scroll_Size - Specifies the vertical size of the draw widget
;	                display in pixels.  Defaults to 512.
;	_EXTRA        - Any allowed Widget_Base keywords may be supplied
;	                through IDL keyword inheritance; these keywords
;	                are passed along to the Widget_Base that provides
;	                the basis of the compound widget.
; COMMENTS:
;	This very simple compound widget simply wraps around a draw widget
;	containing an image to be displayed.  Scrollbars can be used to
;	navigate around the image, allowing a large image to be partially
;	displayed in a small window.
;
;	If the BytScl keyword is specified, the image is loaded into the
;	draw widget graphics window using TVSCL; else TV is used.
;
;	Draw widget events are passed back in a structure CW_IMAGEDRAW,
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
;	Given the id of this compound widget, the procedure CW_IMAGEDRAW_LIMITS
;	will return the array indices of the displayed region.  The syntax
;	is:  CW_ImageDraw_Limits, id, xLL, yLL, xUR, yUR
;
;	As with most compound widgets, it is not recommended that this
;	widget be used as a top level widget.
; MODIFICATION HISTORY:
;	Written by Michael R. Greason, Raytheon ITSS, 02 March 1999.
;	Added a mechanism for returning the indices of the corners of
;	  the displayed region.  MRG, RITSS, 03 March 1999.
;-
; =============================================================================
;
Pro CW_ImageDraw_Limits, id, xLL, yLL, xUR, yUR
;
on_error, 2
;
;			Check arguments.
;
If (n_params() LT 5) Then message, $
	'Syntax: CW_ImageDraw_Limits, id, xLL, yLL, xUR, yUR'
;
;                       Get the state data.
;
stash = Widget_Info(id, /Child)
Widget_Control, stash, Get_UValue=state, /No_Copy
;
;                       Interrogate the draw widget for the current viewport.
;
Widget_Control, state.DrawId, Get_Draw_View=lowleft
;
xLL = lowleft[0]
yLL = lowleft[1]
;
xUR = xLL + state.X_Scroll_Size - 1
yUR = yLL + state.Y_Scroll_Size - 1
;
;                       Reset the state data.
;
Widget_Control, stash, Set_UValue=state, /No_Copy
;
Return
End
;
; =============================================================================
; Paints the draw widget.
;
Pro CW_ImageDraw_DrawImage, state
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
	If (state.Scale) Then tvscl, img $
	                 Else tv,    img
EndIf Else Begin
;
;			Otherwise empty the widget.
;
	erase
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
Function CW_ImageDraw_Value, state, x, y
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
; Centers the displayed portion of the image on a given position.  This will
; not display an out-of-bounds region.
;
Pro CW_ImageDraw_Center, state, xc, yc
;
;			Get the image size.
;
Widget_Control, state.DrawId, Get_Uvalue=img, /No_Copy
sz = size(img)
If (sz[0] EQ 2) Then Begin
;
;			Determine the position of the lower left
;			corner of the displayed area.
;
	xwin = state.X_Scroll_Size
	ywin = state.Y_Scroll_Size
	xsiz = sz[1]
	ysiz = sz[2]
;
	xmax = (xsiz - xwin) > 0
	ymax = (ysiz - ywin) > 0
	xl = ((xc - (xwin / 2)) > 0) < xmax
	yl = ((yc - (ywin / 2)) > 0) < ymax
;
;			Reposition the viewport.
;
	Widget_Control, state.DrawId, Set_Draw_View=[xl, yl]
;
;			Restore the image to the draw widget.
;
EndIf
Widget_Control, state.DrawId, Set_Uvalue=img, /No_Copy
;
Return
End
;
; =============================================================================
; The Widget Realization procedure.  This is called once when the widget
; is realized.  At this time the draw widget has a graphics window assigned
; to it and can be (and is) filled.
;
Pro CW_ImageDraw_Notify, id
;
;                       Get the state data.
;
stash = Widget_Info(id, /Child)
Widget_Control, stash, Get_UValue=state, /No_Copy
;
;                       Paint the draw widget.
;
CW_ImageDraw_DrawImage, state
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
Function CW_ImageDraw_GetValue, id
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
; Assigns the image displayed in the compound widget.
;
Pro CW_ImageDraw_SetValue, id, val
;
;                       Get the state data.
;
stash = Widget_Info(id, /Child)
Widget_Control, stash, Get_UValue=state, /No_Copy
;
;                       Set the value.
;
sz = size(val)
If (sz[0] NE 2) Then img = 0B $
                Else img = val
Widget_Control, state.DrawId, Set_UValue=img, /No_Copy
CW_ImageDraw_DrawImage, state
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
Function CW_ImageDraw_Event, event
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
		typ = ((event.type EQ 0) $
		    OR (event.type EQ 1))
		If (state.BEvent AND typ) Then	$
		    ret = {CW_IMAGEDRAW,	$
			Id:      id,		$
			Top:     event.Top,	$
			Handler: 0L,		$
			Type:    event.Type,	$
			X:       event.X,	$
			Y:       event.Y,	$
			Press:   event.Press,	$
			Release: event.Release,	$
			Clicks:  event.Clicks,	$
			Value:   double(CW_ImageDraw_Value(state, event.X, event.Y))}
;
;				Reposition the view on the mouse
;				click position.
;
		typ = (event.type EQ 0)
		If (state.BCenter and typ) Then $
			CW_ImageDraw_Center, state, event.X, event.Y
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
Function CW_ImageDraw, Img, Parent, BytScl=bscl, 		$
	Button_Events=bevt, Button_Center=bcent,		$
	Retain=ret, X_Scroll_Size=xscrsz, Y_Scroll_Size=yscrsz,	$
	_EXTRA=ext
;
on_error, 2
;
;			Check arguments.
;
If (n_elements(xscrsz) LE 0) Then xscrsz = 1024
If (n_elements(yscrsz) LE 0) Then yscrsz =  512
If (n_elements(ret)    LE 0) Then ret    =    2
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
evfn  = 'CW_ImageDraw_Event'
guvfn = 'CW_ImageDraw_GetValue'
suvpr = 'CW_ImageDraw_SetValue'
nrpr  = 'CW_ImageDraw_Notify'
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
dum = Widget_Base(base)		; Dummy widget used to allow the draw
				; widget user value to be used.
;
drw = Widget_Draw(base,		$
	Retain=ret,		$
	/Scroll,		$
	X_Scroll_Size=xscrsz,	$
	Y_Scroll_Size=yscrsz,	$
	Button_Events=(keyword_set(bevt) OR (keyword_set(bcent))))
;
;			Define the compound widget internal state.
;
state = {DrawId:  drw,			$
         Scale:   keyword_set(bscl),	$
         BEvent:  keyword_set(bevt),	$
	 BCenter: keyword_set(bcent),	$
	 X_Scroll_Size: xscrsz,		$
	 Y_Scroll_Size: yscrsz}
Widget_Control, Widget_Info(base, /Child), Set_UValue=state, /No_Copy
Widget_Control, drw, Set_UValue=val, /No_Copy
;
Return, base
End
