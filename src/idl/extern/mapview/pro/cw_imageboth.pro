;+
; NAME:
;	CW_ImageBoth
; PURPOSE:
;	Defines a compound widget containing a draw widget to be
;	used to display entire images at reduced resolution.
; CALLING SEQUENCE:
;	cw = CW_ImageBoth(Img [, Parent])
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
;	ShowPos       - If present and nonzero, clicking on the reduced
;	                resolution draw widget will cause a box will be
;	                drawn in it showing the portion of the image 
;	                displayed in the other widget.
;	X_Size        - Specifies the horizontal size of the draw widget
;	                displays in pixels.  Defaults to 512.
;	Y_Size        - Specifies the vertical size of the draw widget
;	                displays in pixels.  Defaults to 512.
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
;	Draw widget events are passed back in a structure CW_IMAGEBOTH,
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
;	As with most compound widgets, it is not recommended that this
;	widget be used as a top level widget.
; MODIFICATION HISTORY:
;	Written by Michael R. Greason, Raytheon ITSS, 03 March 1999.
;-
; =============================================================================
; Returns the image displayed in the compound widget.
;
Function CW_ImageBoth_GetValue, id
;
;                       Get the state data.
;
stash = Widget_Info(id, /Child)
Widget_Control, stash, Get_UValue=state, /No_Copy
;
;                       Get the value.
;
Widget_Control, state.ScrollId, Get_Uvalue=val
Widget_Control, state.FullId, Get_Uvalue=val
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
Pro CW_ImageBoth_SetValue, id, val
;
;                       Get the state data.
;
stash = Widget_Info(id, /Child)
Widget_Control, stash, Get_UValue=state, /No_Copy
;
;                       Set the value.
;
sz = size(val)
If (sz[0] NE 2) Then val = 0B
Widget_Control, state.ScrollId, Set_Value=val
Widget_Control, state.FullId,   Set_Value=val
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
Function CW_ImageBoth_Event, event
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
;			Full resolution draw widget events.
;
	'CW_IMAGEDRAW' : Begin
;
;				Return a mouse position event to the
;				parent or application.
;
		If ((event.type EQ 0) OR (event.type EQ 1)) Then $
		    ret = {CW_IMAGEBOTH,	$
			Id:      id,		$
			Top:     event.Top,	$
			Handler: 0L,		$
			Type:    event.Type,	$
			X:       event.X,	$
			Y:       event.Y,	$
			Press:   event.Press,	$
			Release: event.Release,	$
			Clicks:  event.Clicks,	$
			Value:   event.Value}
	End
;
;			Rescaled draw widget events.
;
	'CW_IMAGEFULL' : Begin
;
;				Return a mouse position event to the
;				parent or application.
;
		If ((event.type EQ 0) OR (event.type EQ 1)) Then Begin
		    ret = {CW_IMAGEBOTH,	$
			Id:      id,		$
			Top:     event.Top,	$
			Handler: 0L,		$
			Type:    event.Type,	$
			X:       event.X,	$
			Y:       event.Y,	$
			Press:   event.Press,	$
			Release: event.Release,	$
			Clicks:  event.Clicks,	$
			Value:   event.Value}
;
;				Draw a box in the compressed resolution
;				widget showing where we are.
;
		If ((event.type EQ 0) AND state.ShowPos) Then Begin
		    CW_ImageDraw_Limits, state.ScrollId, xLL,yLL,xUR,yUR
		    CW_ImageFull_Box, state.FullId, xLL,yLL,xUR,yUR
		EndIf
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
Function CW_ImageBoth, Img, Parent, BytScl=bscl, 	$
	Button_Events=bevt, Button_Center=bcent, 	$
	Retain=ret, X_Size=xscrsz, Y_Size=yscrsz,	$
	ShowPos=sp, _EXTRA=ext
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
evfn  = 'CW_ImageBoth_Event'
guvfn = 'CW_ImageBoth_GetValue'
suvpr = 'CW_ImageBoth_SetValue'
;
If (n_params() LT 2) Then Begin
	base = Widget_Base(/Row,		$
			Event_Func=evfn,	$
			Func_Get_Value=guvfn,	$
			Pro_Set_Value=suvpr,	$
			_EXTRA=ext)
EndIf Else Begin
	base = Widget_Base(Parent, /Row,	$
			Event_Func=evfn,	$
			Func_Get_Value=guvfn,	$
			Pro_Set_Value=suvpr,	$
			_EXTRA=ext)
EndElse
;
dum = Widget_Base(base)		; Dummy widget used to allow the draw
				; widget user value to be used.
;
sid = CW_ImageDraw(Img, base,			$
	BytScl=keyword_set(bscl),		$
	Retain=ret,				$
	X_Scroll_Size=xscrsz,			$
	Y_Scroll_Size=yscrsz,			$
	Button_Center=keyword_set(bcent),	$
	Button_Events=keyword_set(bevt))
;
fid = CW_ImageFull(Img, base,			$
	BytScl=keyword_set(bscl),		$
	Retain=ret,				$
	ShowBox=keyword_set(sp),		$
	X_Size=xscrsz,				$
	Y_Size=yscrsz,				$
	Button_Events=keyword_set(bevt))
;
;			Define the compound widget internal state.
;
state = {ScrollId: sid,	$
	 FullId:   fid, $
	 ShowPos:  keyword_set(sp)}
Widget_Control, Widget_Info(base, /Child), Set_UValue=state, /No_Copy
;
Return, base
End
