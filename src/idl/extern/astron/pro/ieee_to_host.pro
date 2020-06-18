pro ieee_to_host, data, IDLTYPE = idltype
;+
; NAME:
;     IEEE_TO_HOST
; PURPOSE:
;     Translate an IDL variable from IEEE-754 to host representation 
; EXPLANATION:
;     The variable is translated from IEEE-754 (as used, for
;     example, in FITS data ), into the host machine architecture.
;
; CALLING SEQUENCE:
;     IEEE_TO_HOST, data, [ IDLTYPE = , ]
;
; INPUT-OUTPUT PARAMETERS:
;     data - any IDL variable, scalar or vector.   It will be modified by
;             IEEE_TO_HOST to convert from IEEE to host representation.  Byte 
;             and string variables are returned by IEEE_TO_HOST unchanged
;
; OPTIONAL KEYWORD INPUTS:
;     IDLTYPE - scalar integer (1-15) specifying the IDL datatype according
;               to the code given by the SIZE function.     This keyword
;               is usually when DATA is a byte array to be interpreted as
;               another datatype (e.g. FLOAT).
;
; EXAMPLE:
;       A 2880 byte array (named FITARR) from a FITS record is to be 
;       interpreted as floating and converted to the host representaton:
;
;       IDL> IEEE_TO_HOST, fitarr, IDLTYPE = 4     
;
; METHOD:
;       The BYTEORDER procedure is called with the appropriate keyword
;
; PROCEDURE CALLS:
;       WHERE_NEGZERO() - Called under VMS prior to V5.1 to check for 
;                         IEEE -0.0 values
;
; MODIFICATION HISTORY:
;      Written, W. Landsman   Hughes/STX   May, 1992
;      Converted to IDL V5.0   W. Landsman   September 1997
;      Under VMS check for IEEE -0.0 values   January 1998
;      VMS now handle -0.0 values under IDL V5.1    July 1998
;      Added new integer datatypes  C. Markwardt/W. Landsman  July 2000
;     
;-
 On_error,2 

 if N_params() EQ 0 then begin
    print,'Syntax - IEEE_TO_HOST, data, [ IDLTYPE = ]'
    return
 endif  

 npts = N_elements( data )
 if npts EQ 0 then $
     message,'ERROR - IDL data variable (first parameter) not defined'

 sz = size(data)
 if not keyword_set( idltype) then idltype = sz[ sz[0]+1]
 
 ; Check for IEEE negative zero values which VMS couldn't handle prior to V5.1

 chk_negzero = (!VERSION.OS EQ 'vms') and (!VERSION.RELEASE LT '5.1')
 if chk_negzero then chk_negzero =  ((idltype GE 4) or (idltype EQ 5) $
                or (idltype EQ 6) or (idltype EQ 9))
 if chk_negzero then wn = WHERE_NEGZERO( data, NCount, /QUIET)


 case idltype of

      1: return                             ;byte

      2: byteorder, data, /NTOHS            ;integer

      3: byteorder, data, /NTOHL            ;long

      4: byteorder, data, /XDRTOF                              ;float

      5: byteorder, data, /XDRTOD                              ;double

      6: byteorder, data, /XDRTOF

      7: return                             ;string

       8: BEGIN				    ;structure

	Ntag = N_tags( data )

	for t=0,Ntag-1 do  begin
          temp = data.(t)
          ieee_to_host, temp
          data.(t) = temp
        endfor 

       END

	9: byteorder, data,/XDRTOD                              ;complex

       12: byteorder, data, /NTOHS

       13: byteorder, data, /NTOHL

       14: if (long(['01'xb,'02'xb,'03'xb,'04'xb],0,1))(0) NE '01020304'x then $
            byteorder, data, /L64swap

       15: if (long(['01'xb,'02'xb,'03'xb,'04'xb],0,1))(0) NE '01020304'x then $  
             byteorder, data, /L64swap

       ELSE: message,'Unrecognized datatype of ' + strtrim(idltype,2)
 ENDCASE

 if chk_negzero then if NCount GT 0 then data[wn] = 0

 return
 end 
