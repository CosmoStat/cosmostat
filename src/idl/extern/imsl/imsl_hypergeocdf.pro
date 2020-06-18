; $Id: //depot/idl/IDL_70/idl_src/libs/imsl/imsl_6.0/lib/imsl_hypergeocdf.pro#1 $
;
; Copyright (c) 1970-2006, VISUAL NUMERICS Inc. All Rights Reserved.
; This software is confidential information which is proprietary to and a
; trade secret of Visual Numerics Inc.  Use, duplication or disclosure is
; subject to the terms of an appropriate license agreement.
;

FUNCTION imsl_hypergeocdf, k, $          ;INPUT Scalar: LONG
              n, $                    ;INPUT Scalar: LONG
              m, $                    ;INPUT Scalar: LONG
              l, $                    ;INPUT Scalar: LONG
              double=double           ;INPUT Scalar ON/OFF flag

@imsl_init.pro   
   ON_ERROR, on_err_action
   ;
   ; Error checking.
   ; Make sure the input argument is a simple array.
   ;
   nargs = n_params()
   IF (nargs EQ 4) THEN BEGIN
      IF (N_ELEMENTS(k) EQ FALSE) THEN message, "K is undefined"
      IF (N_ELEMENTS(n) EQ FALSE) THEN message, "N is undefined"
      IF (N_ELEMENTS(m) EQ FALSE) THEN message, "M is undefined"
      IF (N_ELEMENTS(l) EQ FALSE) THEN message, "L is undefined"
   END ELSE message, "Incorrect number of arguments."
   ;
   ; Setup the parameters for the call to the system function.
   ; 
   ; Input LONG argument(s)
   ;
   k_cvt = IMSL_LONG(k(0))
   n_cvt = IMSL_LONG(n(0))
   m_cvt = IMSL_LONG(m(0))
   l_cvt = IMSL_LONG(l(0))
   IF (KEYWORD_SET(double)) THEN begin
      type = TYP_DOUBLE
      result = double(0.0)
   END ELSE BEGIN
      type = TYP_FLOAT
      result = float(0.0)
   END
   ;
   ; Call the system function.
   ;
   err_status = 0L
   MATHSTAT_149, type, err_status, k_cvt, n_cvt, m_cvt, l_cvt, $
            result
   ;
   ; Return.
   ;
   RETURN, result
END
