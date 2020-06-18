;+
; NAME: 
;    FSTRING
; PURPOSE:
;    Shell around STRING function to fix 1024 size limit on formatting strings
;
; CALLING SEQUENCE:
;    new = fstring(old, [ format, FORMAT = )
;
; INPUTS:
;    OLD = string or number to format, scalar, vector or array
;
; OPTIONAL STRING:
;    FORMAT = scalar string giving format to pass to the STRING() function
;             See restrictions on possible formats below.
; OPTIONAL KEYWORD INPUT:
;    FORMAT  = Format string can alternatively be called as keyword
;
; OUTPUT:
;    FSTRING will return a string with the same dimensions 
;
; RESTRICTIONS:
;    Because FSTRING breaks up the formatting into 1024 element chunks, problems
;    can arise if the number of formatting elements does not evenly divide
;    into 1024.    For example, if format = '(i6,f6.2,e12.6)', (i.e. three
;    formatting elements)  then both the 1023rd and 1024th element will be 
;    formatted as I6.
; EXAMPLE:
;    Create a string array of 10000 uniform random numbers formatted as F6.2
;
;    IDL> a = fstring( randomu(seed,10000), '(f6.2)') 
; REVISION HISTORY:
;     Written W. Landsman (based on program by D. Zarro)  February 2000
;-

function fstring,input,format,format=key_format

 if N_elements(key_format) EQ 1 then form = key_format else $
    if N_elements(format) EQ 1 then  form  = format else form = ''
 if strtrim(form,2) EQ '' then return, string(input)

 np = N_elements(input)
 if np LE 1024 then return,string(input,form=form)
 new = make_array(/string,size = size(input))

; Now format the data in 1024 element chunks

 for i=0L, np-1, 1024 do begin
    i1 = (i + 1023L) < (np-1)
    new[i] = string( input[i:i1], form= form)
 endfor

 return,new 
 end
