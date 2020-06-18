PRO XDUMP, FILENAME=filename, WIN=win, LANDSCAPE=landscape, $
              SCALE=scale, TITLE=title, GIF=gif
;+ 
; NAME: 
;   XDUMP
;
; PURPOSE:
;       save the selected window in a poscript file. The file
;       name is 'idl.ps'. 
;
; CALLING SEQUENCE: 
;   XDUMP, FILENAME=filename, WIN=win, LANDSCAPE=landscape,
;             TITLE=title, SCALE=scale, GIF=gif
;
; OPTIONAL INPUT PARAMETERS: 
;   filename     -- filename of output postscript file (default is idl.ps)
;   win          -- identifier of IDL window (default is active window)
;   landscape    -- set landscape orientation
;   title        -- add the title at the upper left corner of window
;   scale        -- specify a scale applied to the entire graph
;   gif          -- store the data in gif format instead of ps format
;                   Default title in this case idl.gif
;
; EXAMPLES:
;
;   XDUMP
;     Hardcopy the window IDL 0, output file is color.ps
;
;   XDUMP, /LANDSCAPE, SCALE_FACTOR=0.5
;     Same as above but with landscape orientation and a scaling factor of 0.5
;
;   XDUMP, TITLE='my_object'
;     Same first example above but displays also a title  
;
;   XDUMP, FILENAME='my_postscript.ps', WIN=1
;     Hardcopy the window IDL 1, output file is my_postscript.ps
;
; MODIFICATION HISTORY: 
;    11-Mar-1996 A. Claret written.
;    14-Mar-1996 Added title, copyright, scale keywords
;
;-

;------------------------------------------------------------
; parameters check
;------------------------------------------------------------
 
 IF N_PARAMS() LT 0 THEN BEGIN
   PRINT, 'CALLING SEQUENCE: ', $ 
    'XDUMP, filename=filename, win=win, landscape=landscape, ', $
    'title=title, scale=scale, gif=gif'
   GOTO, CLOSING
 ENDIF

;------------------------------------------------------------
; function body
;------------------------------------------------------------

IF keyword_set(win) THEN active_window = win ELSE active_window = !d.window
IF NOT keyword_set(title) THEN title=''

wset, active_window
IF keyword_set(copyright) THEN BEGIN
   xyouts, 0.05, 0.05, /NORMAL, '!3ISOCAM Interactive Analysis', CHARSIZE=1.25
   xyouts, 0.80, 0.05, /NORMAL, STRING("251B)+'!3ESA/CEA-Saclay', CHARSIZE=1.25
ENDIF
xyouts, 0.05, 0.95, /NORMAL, '!6'+title, CHARSIZE=2.0 

a = tvrd(active_window)
IF keyword_set(gif) THEN BEGIN
   IF NOT keyword_set(filename) THEN filename = 'idl.gif'
   write_gif, filename, a
ENDIF ELSE BEGIN
   IF NOT keyword_set(filename) THEN filename = 'idl.ps'
   set_plot, 'PS'
   IF NOT keyword_set(scale) THEN scale = 1.0
   IF keyword_set(landscape) $
      THEN device, file=filename, /color, bits=8, scale_factor=scale, /landscape $
      ELSE device, file=filename, /color, bits=8, scale_factor=scale
   tv, a
   device, /close
   set_plot, 'X'
ENDELSE

;------------------------------------------------------------
; closing
;------------------------------------------------------------
 
 CLOSING:
 
  RETURN
 
 END



