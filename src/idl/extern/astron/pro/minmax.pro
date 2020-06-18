function minmax,array,NAN=nan
;+
; NAME:
;      MINMAX
; PURPOSE:
;      Return a 2 element array giving the minimum and maximum of an array
; EXPLANATION:
;      Using MINMAX() is faster than doing a separate MAX and MIN.
;
; CALLING SEQUENCE:
;      value = minmax( array )
; INPUTS:
;      array - an IDL numeric scalar, vector or array.
;
; OUTPUTS:
;      value = a two element vector, 
;            value[0] = minimum value of array
;            value[1] = maximum value of array
;
; OPTIONAL INPUT KEYWORD:
;      /NAN   - Set this keyword to cause the routine to check for occurrences
;            of the IEEE floating-point value NaN in the input data.  Elements 
;            with the value NaN are treated as missing data.
;
; EXAMPLE:
;      Print the minimum and maximum of an image array, im
; 
;            IDL> print, minmax( im )
;
; PROCEDURE:
;      The MIN function is used with the MAX keyword
;
; REVISION HISTORY:
;      Written W. Landsman                January, 1990
;      Converted to IDL V5.0   W. Landsman   September 1997
;      Added NaN keyword.      M. Buie       June 1998
;-
 On_error,2
 amin = min( array, MAX = amax, NAN=nan)
 return, [ amin, amax ]
 end
