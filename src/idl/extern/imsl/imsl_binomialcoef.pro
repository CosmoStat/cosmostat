; $Id: //depot/idl/IDL_70/idl_src/libs/imsl/imsl_6.0/lib/imsl_binomialcoef.pro#1 $
;
; Copyright (c) 1970-2006, VISUAL NUMERICS Inc. All Rights Reserved.
; This software is confidential information which is proprietary to and a
; trade secret of Visual Numerics Inc.  Use, duplication or disclosure is
; subject to the terms of an appropriate license agreement.
;
FUNCTION imsl_binomialcoef, n, $    ;INPUT Scalar LONG
                   m, $             ;INPUT Scalar LONGt
                   Double=DOUBLE    ;INPUT Scalar ON/OFF flag

@imsl_init.pro
   ON_ERROR, on_err_action
   ;
   ; Error checking.
   ; Make sure the input arguments are a simple arrays.
   ;
   nargs = n_params()
   IF (nargs NE 2) THEN $
         message, "Incorrect number of arguments."
   ;
   ; Decide on what precision to use.
   ;
   type = TYP_FLOAT
   IF (KEYWORD_SET(double) EQ true) THEN type = TYP_DOUBLE
   ;
   ; Setup the parameters for the call to the system function.
   ;
   ; Input LONG argument(s)
   n_cvt = IMSL_LONG(n(0))
   m_cvt = IMSL_LONG(m(0))
   IF (type EQ TYP_DOUBLE) THEN BEGIN
      result = double(0.0)
   END ELSE BEGIN
      result = float(0.0)
   END
   ;
   ; Call the system function.
   ;
   err_status = 0L
   MATHSTAT_250, type, err_status, n_cvt, m_cvt,  $
                 result
   ;
   ; Return.
   ;
   RETURN, result
END
