; $Id: //depot/idl/IDL_70/idl_src/libs/imsl/imsl_6.0/lib/imsl_poissoncdf.pro#1 $
;
; Copyright (c) 1970-2006, VISUAL NUMERICS Inc. All Rights Reserved.
; This software is confidential information which is proprietary to and a
; trade secret of Visual Numerics Inc.  Use, duplication or disclosure is
; subject to the terms of an appropriate license agreement.
;
PRO imsl_Cvt_k_theta, type, k, theta, k_cvt, theta_cvt, result

@imsl_init.pro
   ON_ERROR, on_err_action
   n_elem_k = IMSL_N_ELEMENTS(k)
   n_elem_theta = IMSL_N_ELEMENTS(theta)
   ; CASE 1: One of the arguments is a scalar.
   ;         In this case the result will be a vector whose length
   ;         is determined by the argument with the largest length.
   ;         The special function will be called with pairs made up
   ;         of the scalar and each element of the other argument.
   ;         This case also covers the case when both arguments are
   ;         scalars.
   ;
   IF ((n_elem_k EQ 1) OR (n_elem_theta EQ 1)) THEN BEGIN
      n_elem_result = (n_elem_k > n_elem_theta)
      result = make_array(n_elem_result, type = type)
      IF (n_elem_k EQ 1) THEN BEGIN
         k_cvt = make_array(n_elem_result, type = TYP_MEMINT, value = k(0))
      END ELSE begin
         k_cvt = IMSL_LONG(k)
      END
      IF (n_elem_theta EQ 1) THEN BEGIN
         theta_cvt = make_array(n_elem_result, type = type, value = theta(0))
      END ELSE begin
         IF (type EQ TYP_DOUBLE) THEN theta_cvt = double(theta) $
           ELSE theta_cvt = float(theta)
      END
   END ELSE BEGIN
   ; CASE 2: Both arguments are vectors
   ;         In this case the result will be a vector whose length
   ;         is determined by the argument with the shortest length.
   ;         The special function will be called with pairs made up
   ;         the i-th element of the each vector.
   ;
      n_elem_result = (n_elem_k < n_elem_theta)
      result = make_array(n_elem_result, type = type)
      IF (n_elem_theta LT n_elem_k) THEN BEGIN
         k_cvt = make_array(n_elem_result, type = TYP_MEMINT)
         k_cvt(*) = k(0:n_elem_result-1)
         IF (type EQ TYP_DOUBLE) THEN theta_cvt = double(theta) $
           ELSE theta_cvt = float(theta)
      END ELSE BEGIN
         theta_cvt = make_array(n_elem_result, type = type)
         theta_cvt(*) = theta(0:n_elem_result-1)
         k_cvt = IMSL_LONG(k)
      END
   END
   IF (n_elem_result EQ 1) THEN result = result(0)
END

FUNCTION imsl_poissoncdf, k, $          ;INPUT Scalar or 1-D array: LONG
              theta, $                ;INPUT scalar or 1-D array: floating point
              double=double           ;INPUT Scalar ON/OFF flag
@imsl_init.pro
   ON_ERROR, on_err_action
   ;
   ; Error checking.
   ; Make sure the input argument is a simple array.
   ;
   nargs = n_params()
   IF (nargs EQ 2) THEN BEGIN
      IF (N_ELEMENTS(k) EQ FALSE) THEN message, "K is undefined"
      IF (N_ELEMENTS(theta) EQ FALSE) THEN message, "THETA is undefined"
      size_k = IMSL_SIZE(k)
      size_theta = IMSL_SIZE(theta)
      IF (size_k(0) GT 1) THEN $
        message, "K must be either a scalar or a 1-D array."
      IF (size_theta(0) GT 1) THEN $
        message, "THETA must be either a scalar or a 1-D array."
   END ELSE message, "Incorrect number of arguments."
   ;
   ; Decide on what precision to use.
   ; Use the highest precision of the input argument(s).
   ;
   type = TYP_FLOAT
   IF (size_theta(N_ELEMENTS(size_theta)-2) GT TYP_FLOAT) THEN type = TYP_DOUBLE
   IF (KEYWORD_SET(double) EQ true) THEN type = TYP_DOUBLE
   ;
   ; Setup the parameters for the call to the system function.
   ;
   ;
   ; Input LONG argument(s)
   ;
   ; Floating point arguments and keywords
   imsl_Cvt_k_theta, type, k, theta, k_cvt, theta_cvt, result
   ;
   ; Call the system function.
   ;
   err_status = 0L
   MATHSTAT_172, type, err_status, k_cvt, theta_cvt, $
            IMSL_N_ELEMENTS(result), $
            result
   ;
   ; Return.
   ;
   RETURN, result
END
