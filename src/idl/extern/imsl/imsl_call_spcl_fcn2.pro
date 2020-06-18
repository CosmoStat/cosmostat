; $Id: //depot/idl/IDL_70/idl_src/libs/imsl/imsl_6.0/lib/imsl_call_spcl_fcn2.pro#1 $
;
; Copyright (c) 1970-2006, VISUAL NUMERICS Inc. All Rights Reserved.
; This software is confidential information which is proprietary to and a
; trade secret of Visual Numerics Inc.  Use, duplication or disclosure is
; subject to the terms of an appropriate license agreement.
;
PRO imsl_Cvt_x_y, type, x, y, x_cvt, y_cvt, result

@imsl_init.pro
   ON_ERROR, on_err_action
   ; Initialization of local variables.
   n_elem_x = IMSL_N_ELEMENTS(x)
   n_elem_y = IMSL_N_ELEMENTS(y)
   ; CASE 1: One of the arguments is a scalar.
   ;         In this case the result will be a vector whose length
   ;         is determined by the argument with the largest length.
   ;         The special function will be called with pairs made up
   ;         of the scalar and each element of the other argument.
   ;         This case also covers the case when both arguments are
   ;         scalars.
   ;
   IF ((n_elem_x EQ 1) OR (n_elem_y EQ 1)) THEN BEGIN
      n_elem_result = (n_elem_x > n_elem_y)
      result = make_array(n_elem_result, type = type)
      IF (n_elem_x EQ 1) THEN BEGIN
         x_cvt = make_array(n_elem_result, type = type, value = x(0))
      END ELSE begin
         IF (type EQ TYP_DOUBLE) THEN x_cvt = double(x) ELSE x_cvt = float(x)
      END
      IF (n_elem_y EQ 1) THEN BEGIN
         y_cvt = make_array(n_elem_result, type = type, value = y(0))
      END ELSE begin
         IF (type EQ TYP_DOUBLE) THEN y_cvt = double(y) ELSE y_cvt = float(y)
      END
   END ELSE BEGIN
   ; CASE 2: Both arguments are vectors
   ;         In this case the result will be a vector whose length
   ;         is determined by the argument with the shortest length.
   ;         The special function will be called with pairs made up
   ;         the i-th element of the each vector.
   ;
      n_elem_result = (n_elem_x < n_elem_y)
      result = make_array(n_elem_result, type = type)
      IF (n_elem_y LT n_elem_x) THEN BEGIN
         x_cvt = make_array(n_elem_result, type = type)
         x_cvt(*) = x(0:n_elem_result-1)
         IF (type EQ TYP_DOUBLE) THEN y_cvt = double(y) ELSE y_cvt = float(y)
      END ELSE BEGIN
         y_cvt = make_array(n_elem_result, type = type)
         y_cvt(*) = y(0:n_elem_result-1)
         IF (type EQ TYP_DOUBLE) THEN x_cvt = double(x) ELSE x_cvt = float(x)
      END
   END
   IF (n_elem_result EQ 1) THEN result = result(0)
END

FUNCTION imsl_call_spcl_fcn2, x, $         ;INPUT 1-D array: floating point
              y, $                    ;INPUT 1-D array: floating point
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
   IF (nargs EQ 2) THEN BEGIN
      IF (N_ELEMENTS(x) EQ FALSE) THEN message, "X is undefined"
      IF (N_ELEMENTS(y) EQ FALSE) THEN message, "Y is undefined"
      size_x = IMSL_SIZE(x)
      size_y = IMSL_SIZE(y)
      IF (size_x(0) GT 1) THEN $
        message, "X must be either a scalar or a 1-D array."
      IF (size_y(0) GT 1) THEN $
        message, "Y must be either a scalar or a 1-D array."
   END ELSE message, "Incorrect number of arguments."
   ;
   ; Decide on what precision to use.
   ; Use the highest precision of the input argument(s).
   ;
   type = TYP_FLOAT
   size_x = IMSL_SIZE(x)
   size_y = IMSL_SIZE(y)
   IF (size_x(N_ELEMENTS(size_x)-2) GT TYP_FLOAT) THEN type = TYP_DOUBLE
   IF (size_y(N_ELEMENTS(size_y)-2) GT TYP_FLOAT) THEN type = TYP_DOUBLE
   IF (KEYWORD_SET(double) EQ true) THEN type = TYP_DOUBLE
   ;
   ; Setup the parameters for the call to the system function.
   ;
   ;
   ; Input LONG argument(s)
   nx = IMSL_N_ELEMENTS(x)
   ny = IMSL_N_ELEMENTS(y)
   ;
   ; Floating point arguments and keywords
   imsl_Cvt_x_y, type, x, y, x_cvt, y_cvt, result
   ;
   ; Call the system function.
   ;
   err_status = 0L

   MATHSTAT_135, type, err_status, x_cvt, $
            IMSL_N_ELEMENTS(result), $
            IMSL_LONG(KEYWORD_SET(inverse)), $
            fcn_name, result, y_cvt
   ;
   ; Return.
   ;
   RETURN, result
END
