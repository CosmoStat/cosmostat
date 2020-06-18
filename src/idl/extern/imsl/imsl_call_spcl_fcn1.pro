; $Id: //depot/idl/IDL_70/idl_src/libs/imsl/imsl_6.0/lib/imsl_call_spcl_fcn1.pro#1 $
;
; Copyright (c) 1970-2006, VISUAL NUMERICS Inc. All Rights Reserved.
; This software is confidential information which is proprietary to and a
; trade secret of Visual Numerics Inc.  Use, duplication or disclosure is
; subject to the terms of an appropriate license agreement.
;

FUNCTION imsl_call_spcl_fcn1, x, $         ;INPUT 1-D array: floating point
              fcn_name=fcn_name, $    ;INPUT Scalar STRING
              double=double, $        ;INPUT Scalar ON/OFF flag
              inverse=inverse         ;INPUT Scalar ON/OFF flag

@imsl_init.pro
   ON_ERROR, on_err_action
   ;
   ; Error checking.
   ; Make sure the input argument is a simple array.
   ;
   nargs = n_params()
   IF (nargs EQ 1) THEN BEGIN
      IF (N_ELEMENTS(x) EQ FALSE) THEN message, "X is undefined"
      size_x = IMSL_SIZE(x)
      IF (size_x(0) GT 1) THEN $
        message, "X must be either a scalar or a 1-D array."
   END ELSE message, "Incorrect number of arguments."
   ;
   ; Decide on what precision to use.
   ; Use the highest precision of the input argument(s).
   ;
   type = TYP_FLOAT
   IF (size_x(N_ELEMENTS(size_x)-2) GT TYP_FLOAT) THEN type = TYP_DOUBLE
   IF (KEYWORD_SET(double) EQ true) THEN type = TYP_DOUBLE
   ;
   ; Setup the parameters for the call to the system function.
   ;
   ;
   ; Input LONG argument(s)
   nx = IMSL_N_ELEMENTS(x)
   ;
   ; Floating point arguments and keywords
   IF (type EQ TYP_DOUBLE) THEN BEGIN
      ;
      ; Input
      x_cvt = double(x)
      ; Result
      ;
      IF (N_ELEMENTS(x) EQ 1) THEN result  = double(0.0) ELSE $
        result = make_array(size = IMSL_SIZE(x_cvt))
   END ELSE BEGIN
      ;
      ; Input
      x_cvt = float(x)
      ; Result
      ;
      IF (N_ELEMENTS(x) EQ 1) THEN result  = float(0.0) ELSE $
        result = make_array(size = IMSL_SIZE(x_cvt))
   END
   ;
   ; Call the system function.
   ;
   err_status = 0L
   MATHSTAT_135, type, err_status, x_cvt, nx, IMSL_LONG(KEYWORD_SET(inverse)), $
            fcn_name, result
   ;
   ; Return.
   ;
   RETURN, result
END
