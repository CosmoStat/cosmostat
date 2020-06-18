; $Id: //depot/idl/IDL_70/idl_src/libs/imsl/imsl_6.0/lib/imsl_binomialpdf.pro#1 $
;
; Copyright (c) 1970-2006, VISUAL NUMERICS Inc. All Rights Reserved.
; This software is confidential information which is proprietary to and a
; trade secret of Visual Numerics Inc.  Use, duplication or disclosure is
; subject to the terms of an appropriate license agreement.
;

FUNCTION imsl_binomialpdf, k, $       ;INPUT Scalar: LONG
              n, $                    ;INPUT Scalar: LONG
              p, $                    ;INPUT Scalar: LONG
              double=double           ;INPUT Scalar ON/OFF flag
   
@imsl_init.pro
   ON_ERROR, on_err_action
   ;
   ; Error checking.
   ; Make sure the input argument is a simple array.
   ;
   nargs = n_params()
   IF (nargs EQ 3) THEN BEGIN
      IF (N_ELEMENTS(k) EQ FALSE) THEN message, "K is undefined"
      IF (N_ELEMENTS(n) EQ FALSE) THEN message, "N is undefined"
      IF (N_ELEMENTS(p) EQ FALSE) THEN message, "P is undefined"
   END ELSE message, "Incorrect number of arguments."
   ;
   ; Decide on what precision to use.
   ; Use the highest precision of the input argument(s).
   ;
   type = TYP_FLOAT
   size_p = IMSL_SIZE(p)
   IF (size_p(N_ELEMENTS(size_p)-2) GT TYP_FLOAT) THEN type = TYP_DOUBLE
   IF (KEYWORD_SET(double) EQ true) THEN type = TYP_DOUBLE
   ;
   ; Setup the parameters for the call to the system function.
   ; 
   ; Input LONG argument(s)
   ;
   k_cvt = IMSL_LONG(k(0))
   n_cvt = IMSL_LONG(n(0))
   
   IF (type EQ TYP_DOUBLE) THEN BEGIN
      p_cvt = double(p(0))
      result = double(0.0)
   END ELSE BEGIN
      p_cvt = float(p(0))
      result = float(0.0)
   END
   ;
   ; Call the system function.
   ;
   err_status = 0L
   MATHSTAT_310, type, err_status, k_cvt, n_cvt, p_cvt,  $
            result
   ;
   ; Return.
   ;
   RETURN, result
END
