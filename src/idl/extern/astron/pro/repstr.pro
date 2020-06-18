function repstr,obj,in,out
;+
; NAME:
;	REPSTR
; PURPOSE:
;	Replace all occurences of one substring by another.
; EXPLANATION:
;	Meant to emulate the string substitution capabilities of text editors
; CALLING SEQUENCE:
;	result = repstr( obj, in, out )
;
; INPUT PARAMETERS:
;	obj    = object string for editing, scalar or array
;	in     = substring of 'obj' to be replaced, scalar 
;
; OPTIONAL INPUT PARMETER:
;	out    = what 'in' is replaced with, scalar.   If not supplied
;		then out = '', i.e. 'in' is not replaced by anything. 
;
; OUTPUT PARAMETERS:
;	Result returned as function value.  Input object string
;	not changed unless assignment done in calling program.
;
; PROCEDURE:
;	Searches for 'in', splits 'obj' into 3 pieces, reassembles
;	with 'out' in place of 'in'.  Repeats until all cases done.
;
; EXAMPLE:
;	If a = 'I am what I am' then print,repstr(a,'am','was')
;	will give 'I was what I was'.
;
; MODIFICATION HISTORY:
;	Written by Robert S. Hill, ST Systems Corp., 12 April 1989.
;	Accept vector object strings, W. Landsman   HSTX,   April, 1996
;	Converted to IDL V5.0   W. Landsman   September 1997
;-
 On_error,2
 if N_params() LT 3 then begin
	print,'Syntax - result = REPSTR( obj, in, out )'
	return, obj
 endif

 if N_elements(out) EQ 0 then out = ''
 l1 = strlen(in)
 l2 = strlen(out)
 Nstring = N_elements(obj)
 object = obj
 for i= 0 ,Nstring-1 do begin
 last_pos = 0
 lo = 9999
 pos = 0
 while ( pos LT lo-l1) and (pos GE 0) do begin
   lo = strlen(object[i])
   pos = strpos(object[i],in,last_pos)
   if (pos GE 0) then begin
	      first_part = strmid(object[i],0,pos)
	      last_part  = strmid(object[i],pos+l1,9999)
	      object[i] = first_part + out + last_part
   endif 
  last_pos = pos + l2
 endwhile
 endfor

 return,object

 end
