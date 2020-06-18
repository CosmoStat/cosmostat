;+
; NAME:
;	ScrTV
; PURPOSE:
;	Displays an image in a scrollable draw widget.  This allows the
;	user to examine a full resolution portion of the image, using
;	the scrollbars to navigate around.
; CALLING SEQUENCE:
;	ScrTV, img [, cw]
; INPUTS:
;	img - The image to be initially displayed in the application.
; OPTIONAL OUTPUTS:
;	cw  - The widget id of the compound widget.
; OTPIONAL INPUT KEYWORDS:
;	/Button_Center - If present and nonzero, the view is centered on
;	                to the position of a mouse button click.
;	/Button_Events - If present and nonzero, the position of a mouse
;	                click in the draw view is written to the screen
;	                with its corresponding image value.
;	/BytScl        - If present and nonzero, the image is bytescaled
;	                before it is displayed.
;	/Full          - If present and nonzero, a compound widget that 
;	                displays the entire (compressed) image without
;	                scrollbars is added.
;	Retain        - Allows the user to specify the backing store to
;	                use.  Defaults to 2 -- IDL provides backing store.
;	ShowPos       - If present and nonzero, clicking on the reduced
;	                resolution draw widget will cause a box will be
;	                drawn in it showing the portion of the image 
;	                displayed in the other widget.
;	X_Scroll_Size - Specifies the horizontal size of the draw widget
;	                display in pixels.  Defaults to 1024 or 512.
;	Y_Scroll_Size - Specifies the vertical size of the draw widget
;	                display in pixels.  Defaults to 512.
; COMMENTS:
;	The CW_IMAGEDRAW compound widget provides nearly all the functionality.
;	This application adds the position/value display feature, and a
;	Done button at the bottom of the widget application so that the
;	application can be gracefully terminated.
;
;	This is not a blocking application, so IDL returns to the commandline
;	prompt without waiting for this application to terminate.
;
;	If the BytScl keyword is specified, the image is loaded into the
;	draw widget graphics window using TVSCL; else TV is used.
;
;	If Full is specified, both display widgets are given the same size,
;	and the default X size becomes 512.  The two widgets are displayed
;	side by side using the CW_ImageBoth compound widget.
;
;	If the user returns the compound widget id (the cw argument), the
;	Widget_Control procedure can be used to change the image displayed.
;	Simply change the value of the widget:
;		Widget_Control, cw, Set_Value=newimage
; MODIFICATION HISTORY:
;	Written by Michael R. Greason, Raytheon ITSS, 02 March 1999.
;	Full image/low resolution display added.  MRG, RITSS, 03 March 1999.
;-
; =============================================================================
; The application event handler.
;
Pro ScrTV_Event, event
;
tn  = strupcase(strtrim(Tag_Names(event, /Structure_Name), 2))
;
;		Process the event.
;
Case (tn) Of
;
;			Button events.
;
	'WIDGET_BUTTON' : Begin
		Widget_Control, event.id, Get_UValue=cmnd
;
;				The Done button.
;
		If (cmnd EQ 'Done') Then Widget_Control, event.top, /Destroy
;
	End
;
;			Draw events.
;
	'CW_IMAGEDRAW' : Begin
		If (event.type eq 0) Then $
			print, event.X, event.Y, event.Value
	End
;
	'CW_IMAGEBOTH' : Begin
		If (event.type eq 0) Then $
			print, event.X, event.Y, event.Value
	End
;
	Else : ret = 0
EndCase
;
Return
End
;
; =============================================================================
; The entry point into this procedure/application.
;
Pro ScrTV, img, cw, BytScl=bscl, Button_Events=bevt, Button_Center=bcent, $
	Full=ful, Retain=ret, X_Scroll_Size=xsz, Y_Scroll_Size=ysz, ShowPos=sp
;
on_error, 2
;
;			Check arguments.
;
If (n_params() LT 1) Then message, 'Syntax:  ScrTV, img [, cw]'
;
bev = keyword_set(bevt)
bct = keyword_set(bcent)
If (n_elements(ret) LE 0) Then ret = 2
;
If (keyword_set(ful)) Then xsiz =  512 $
                      Else xsiz = 1024
If (n_elements(xsz) GT 0) Then Begin
	If (xsz GT 0) Then xsiz = xsz
EndIf
;
ysiz = 512
If (n_elements(ysz) GT 0) Then Begin
	If (ysz GT 0) Then ysiz = ysz
EndIf
;
base = Widget_Base(/Column)
;
If (keyword_set(ful)) Then Begin
	cw   = CW_ImageBoth(img, base,		$
		BytScl=keyword_set(bscl),	$
		Button_Center=bct,		$
		Button_Event=bev,		$
		Retain=ret,			$
		ShowPos=keyword_set(sp),	$
		X_Size=xsiz,		 	$
		Y_Size=ysiz)
EndIf Else Begin
	cw   = CW_ImageDraw(img, base,		$
		BytScl=keyword_set(bscl),	$
		Button_Center=bct,		$
		Button_Event=bev,		$
		Retain=ret,			$
		X_Scroll_Size=xsiz,		$
		Y_Scroll_Size=ysiz)
EndElse
;
but  = Widget_Button(base, Value='Done', UValue='Done')
;
Widget_Control, base, /Realize
Xmanager, 'ScrTV', base, /No_Block
;
Return
End
