; $Id: //depot/idl/IDL_70/idl_src/libs/imsl/imsl_6.0/lib/imsl_elrj.pro#1 $
;
; Copyright (c) 1970-2006, VISUAL NUMERICS Inc. All Rights Reserved.
; This software is confidential information which is proprietary to and a
; trade secret of Visual Numerics Inc.  Use, duplication or disclosure is
; subject to the terms of an appropriate license agreement.
;

FUNCTION imsl_elrj, x, $          	      ;INPUT Scalar: LONG
              y, $                    ;INPUT Scalar: LONG
              z, $                    ;INPUT Scalar: LONG
              rho, $                  ;INPUT Scalar: LONG
              double=double           ;INPUT Scalar ON/OFF flag

@imsl_init.pro
   ON_ERROR, on_err_action
   ;
   ; Error checking.
   ; Make sure the input argument is a simple array.
   ;
   nargs = n_params()
   IF (nargs EQ 4) THEN BEGIN
      IF (N_ELEMENTS(x) EQ FALSE) THEN message, "X is undefined"
      IF (N_ELEMENTS(y) EQ FALSE) THEN message, "Y is undefined"
      IF (N_ELEMENTS(z) EQ FALSE) THEN message, "Z is undefined"
      IF (N_ELEMENTS(rho) EQ FALSE) THEN message, "RHO is undefined"
   END ELSE message, "Incorrect number of arguments."
   ;
   ; Decide on what precision to use.
   ; Use the highest precision of the input argument(s).
   ;
   type = TYP_FLOAT
   size_tmp = IMSL_SIZE(x)
   IF (size_tmp(N_ELEMENTS(size_tmp)-2) GT TYP_FLOAT) THEN type = TYP_DOUBLE
   size_tmp = IMSL_SIZE(y)
   IF (size_tmp(N_ELEMENTS(size_tmp)-2) GT TYP_FLOAT) THEN type = TYP_DOUBLE
   size_tmp = IMSL_SIZE(z)
   IF (size_tmp(N_ELEMENTS(size_tmp)-2) GT TYP_FLOAT) THEN type = TYP_DOUBLE
   size_tmp = IMSL_SIZE(rho)
   IF (size_tmp(N_ELEMENTS(size_tmp)-2) GT TYP_FLOAT) THEN type = TYP_DOUBLE
   IF (KEYWORD_SET(double) EQ true) THEN type = TYP_DOUBLE
   ;
   ; Setup the parameters for the call to the system function.
   ; 
   
   IF (type EQ TYP_DOUBLE) THEN BEGIN
      x_cvt = double(x(0))
      y_cvt = double(y(0))
      z_cvt = double(z(0))
      rho_cvt = double(rho(0))
      result = double(0.0)
   END ELSE BEGIN
      x_cvt = float(x(0))
      y_cvt = float(y(0))
      z_cvt = float(z(0))
      rho_cvt = float(rho(0))
      result = float(0.0)
   END
   ;
   ; Call the system function.
   ;
   err_status = 0L
   MATHSTAT_213, type, err_status, x_cvt, y_cvt, z_cvt, rho_cvt,  $
            result
   ;
   ; Return.
   ;
   RETURN, result
END
