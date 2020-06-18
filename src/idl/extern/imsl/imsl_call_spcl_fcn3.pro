; $Id: //depot/idl/IDL_70/idl_src/libs/imsl/imsl_6.0/lib/imsl_call_spcl_fcn3.pro#1 $
;
; Copyright (c) 1970-2006, VISUAL NUMERICS Inc. All Rights Reserved.
; This software is confidential information which is proprietary to and a
; trade secret of Visual Numerics Inc.  Use, duplication or disclosure is
; subject to the terms of an appropriate license agreement.
;
PRO imsl_Cvt_1scalar, arg1, arg2, arg3, arg1_cvt, arg2_cvt, arg3_cvt, $
               n_elem_result, type

@imsl_init.pro
   ON_ERROR, on_err_action

   arg1_cvt = make_array(n_elem_result, type = type, value = arg1(0))
   arg2_cvt = arg2(0:n_elem_result-1)
   arg3_cvt = arg3(0:n_elem_result-1)
   IF (type EQ TYP_DOUBLE) THEN BEGIN
      arg2_cvt = double(arg2(0:n_elem_result-1))
      arg3_cvt = double(arg3(0:n_elem_result-1))
   END ELSE BEGIN
      arg2_cvt = float(arg2(0:n_elem_result-1))
      arg3_cvt = float(arg3(0:n_elem_result-1))
   END
END


PRO imsl_Cvt_2scalars, arg1, arg2, arg3, arg1_cvt, arg2_cvt, arg3_cvt, $
               n_elem_result, type
@imsl_init.pro
   ON_ERROR, on_err_action

   arg1_cvt = make_array(n_elem_result, type = type, value = arg1(0))
   arg2_cvt = make_array(n_elem_result, type = type, value = arg2(0))
   arg3_cvt = arg3(0:n_elem_result-1)
   IF (type EQ TYP_DOUBLE) THEN BEGIN
      arg3_cvt = double(arg3)
   END ELSE BEGIN
      arg3_cvt = float(arg3)
   END
END


PRO imsl_Cvt_x_y_z, type, x, y, z, x_cvt, y_cvt, z_cvt, result

@imsl_init.pro
   ON_ERROR, on_err_action

   n_elem_x = IMSL_N_ELEMENTS(x)
   n_elem_y = IMSL_N_ELEMENTS(y)
   n_elem_z = IMSL_N_ELEMENTS(z)
   lengths = [n_elem_x, n_elem_y, n_elem_z]
   lengths = lengths(sort(lengths)) ;Sort the lengths
   ; CASE 1: Exactly One of the arguments is a scalar.
   ;         The length of the result will be the length of the
   ;         shortest vector argument.
   ; CASE 2: Exactly Two of the arguments are scalars.
   ;         The special function will be called with triples made up
   ;         of the scalars and each element of the other argument.
   ; CASE 3: Exactly Three of the arguments are scalars.
   ;         Convert each scalar and the result is a scalar.
   ; CASE 4: Exactly none (0) of the arguments are scalars.
   ;         The result will be the length of the shortest vector.
   ;
   CASE ((n_elem_x EQ 1) + (n_elem_y EQ 1) + (n_elem_z EQ 1)) OF
      0: BEGIN                  ; Exactly 0 scalars.
         n_elem_result =lengths(0)
         result = MAKE_ARRAY(n_elem_result, type = type)
         x_cvt = x(0:n_elem_result-1)
         y_cvt = y(0:n_elem_result-1)
         z_cvt = z(0:n_elem_result-1)
         IF (type EQ TYP_DOUBLE) THEN BEGIN
            x_cvt = double(x(0:n_elem_result-1))
            y_cvt = double(y(0:n_elem_result-1))
            z_cvt = double(z(0:n_elem_result-1))
         END ELSE BEGIN
            x_cvt = float(x(0:n_elem_result-1))
            y_cvt = float(y(0:n_elem_result-1))
            z_cvt = float(z(0:n_elem_result-1))
         END
      END
      1: BEGIN                  ; Exactly 1 scalar.
         n_elem_result =lengths(1)
         result = make_array(n_elem_result, type = type)
         IF (n_elem_x EQ 1) THEN $
           imsl_Cvt_1scalar, x, y, z, x_cvt, y_cvt, z_cvt, n_elem_result, type
         IF (n_elem_y EQ 1) THEN $
           imsl_Cvt_1scalar, y, x, z, y_cvt, x_cvt, z_cvt, n_elem_result, type
         IF (n_elem_z EQ 1) THEN $
           imsl_Cvt_1scalar, z, y, x, z_cvt, y_cvt, x_cvt, n_elem_result, type
      END
      2: BEGIN                  ; Exactly 2 scalars.
         n_elem_result =lengths(2)
         result = make_array(n_elem_result, type = type)
         IF (n_elem_x NE 1) THEN $
           imsl_Cvt_2scalars, y, z, x, y_cvt, z_cvt, x_cvt, n_elem_result, type
         IF (n_elem_y NE 1) THEN $
           imsl_Cvt_2scalars, x, z, y, x_cvt, z_cvt, y_cvt, n_elem_result, type
         IF (n_elem_z NE 1) THEN $
           imsl_Cvt_2scalars, x, y, z, x_cvt, y_cvt, z_cvt, n_elem_result, type
      END
      3: BEGIN                  ; Exactly 3 scalars.
         n_elem_result = IMSL_1
         IF (type EQ TYP_DOUBLE) THEN BEGIN
            result= double(0.0)
            x_cvt = double(x(0))
            y_cvt = double(y(0))
            z_cvt = double(z(0))
         END ELSE BEGIN
            result = float(0.0)
            x_cvt = float(x(0))
            y_cvt = float(y(0))
            z_cvt = float(z(0))
         END
      END
   END
   IF (n_elem_result EQ 1) THEN result = result(0)
END

FUNCTION imsl_call_spcl_fcn3, x, $         ;INPUT 1-D array: floating point
              y, $                    ;INPUT 1-D array: floating point
              z, $                    ;INPUT 1-D array: floating point
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
   IF (nargs EQ 3) THEN BEGIN
      IF (N_ELEMENTS(x) EQ FALSE) THEN message, "X is undefined"
      IF (N_ELEMENTS(y) EQ FALSE) THEN message, "Y is undefined"
      IF (N_ELEMENTS(z) EQ FALSE) THEN message, "Z is undefined"
      size_x = IMSL_SIZE(x)
      size_y = IMSL_SIZE(y)
      size_z = IMSL_SIZE(z)
      IF (size_x(0) GT 1) THEN $
        message, "X must be either a scalar or a 1-D array."
      IF (size_y(0) GT 1) THEN $
        message, "Y must be either a scalar or a 1-D array."
      IF (size_z(0) GT 1) THEN $
        message, "Z must be either a scalar or a 1-D array."
   END ELSE message, "Incorrect number of arguments."
   ;
   ; Decide on what precision to use.
   ; Use the highest precision of the input argument(s).
   ;
   type = TYP_FLOAT
   IF (size_x(N_ELEMENTS(size_x)-2) GT TYP_FLOAT) THEN type = TYP_DOUBLE
   IF (size_y(N_ELEMENTS(size_y)-2) GT TYP_FLOAT) THEN type = TYP_DOUBLE
   IF (size_z(N_ELEMENTS(size_z)-2) GT TYP_FLOAT) THEN type = TYP_DOUBLE
   IF (KEYWORD_SET(double) EQ true) THEN type = TYP_DOUBLE
   ;
   ; Setup the parameters for the call to the system function.
   ;
   ; Floating point arguments and keywords
   imsl_Cvt_x_y_z, type, x, y, z, x_cvt, y_cvt, z_cvt, result
   ;
   ; Call the system function.
   ;
   err_status = 0L
   MATHSTAT_135, type, err_status, x_cvt, $
            IMSL_N_ELEMENTS(result), $
            IMSL_LONG(KEYWORD_SET(inverse)), $
            fcn_name, result, y_cvt, z_cvt
   ;
   ; Return.
   ;
   RETURN, result
END
