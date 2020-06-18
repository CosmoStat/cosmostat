pro prt_precision, eps=eps, double=double

;Aug 2008, Modified by An.R, removed redundant _extra keyword
;Written by Anais Rassat, July 2008.

; ***************// HELP BEGIN //**************
; PURPOSE:  This routine is designed to print to your screen the
;           precision of your machine in either single or double
;           precision.
; INPUT:    None
; OPTIONAL INPUT:
;           double: set this to 1 to obtain the machine's double
;                   precision
; OUTPUT:   to screen: the precision of your machine, detailing single
;                   or double precision
; OPTIONAL OUTPUT:
;           eps: precision
; ---
; EXAMPLE 1:
;          > prt_precision
; IDL returns to screen:
;          Accuracy in single precision is:   1.19209e-07
; 
; EXAMPLE 2: 
;          > prt_precision, /double
; IDL returns to screen:
;          Accuracy in double precision is:    2.2204460e-16         
; ***************// HELP END //**************

  eps = (machar(double=double)).eps
  
if keyword_set(double) then prec = 'double' else prec = 'single'
print, 'Accuracy in '+prec+' precision is: ', eps

end
