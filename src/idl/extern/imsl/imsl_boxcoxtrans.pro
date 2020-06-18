; $Id: //depot/idl/IDL_70/idl_src/libs/imsl/imsl_6.0/lib/imsl_boxcoxtrans.pro#1 $
;
; Copyright (c) 1970-2006, VISUAL NUMERICS Inc. All Rights Reserved.
; This software is confidential information which is proprietary to and a
; trade secret of Visual Numerics Inc.  Use, duplication or disclosure is
; subject to the terms of an appropriate license agreement.
;
FUNCTION imsl_boxcoxtrans, z, $           ;INPUT 1-D array: floating point 
           power,  $                      ;INPUT Scalar floating point
           inverse=inverse, $             ;INPUT Scalar ON/OFF flag 
           s=s, $                         ;INPUT Scalar floating point 
           double=double                  ;INPUT Scalar ON/OFF flag
           
@imsl_init.pro
   ON_ERROR, on_err_action
   ;
   ; Error checking:
   ;  - Z must be either a 1-D array. (set n = length)
   ;          
   nargs = n_params()
   IF (nargs NE 2)  THEN message, 'Incorrect number of arguments.'
     
   size_z = IMSL_LONG(size(z))
   IF (size_z(0) NE 1) THEN message, 'Z must be a 1-D array.'
   n = IMSL_LONG(size_z(1))
   ;
   ; Decide on what precision to use.
   ;
   type = TYP_FLOAT
   IF (size_z(N_ELEMENTS(size_z)-2) EQ  TYP_DOUBLE) THEN type = TYP_DOUBLE
   IF (KEYWORD_SET(double) EQ true) THEN type = TYP_DOUBLE
   ;
   ; Setup the parameters for the call to the system function.
   ;
   IF (KEYWORD_SET(inverse)) THEN inverse_cvt = IMSL_1 ELSE inverse_cvt = IMSL_0
   ;
   ; Floating point arguments and keywords
   ;
   IF (type EQ TYP_DOUBLE) THEN BEGIN
      result = dblarr(n)
      z_cvt = double(z)
      power_cvt = double(power(0))
      IF (KEYWORD_SET(s)) THEN s_cvt  =  DOUBLE(s(0)) ELSE s_cvt  =  DOUBLE(0.0)
   END ELSE BEGIN
      result = fltarr(n)
      z_cvt = float(z)
      power_cvt = float(power(0))
      IF (KEYWORD_SET(s)) THEN s_cvt  =  FLOAT(s(0)) ELSE s_cvt  =  FLOAT(0.0)
   END
   ;
   ; Call the system function.
   ;
   err_status = 0L
   MATHSTAT_251, type, err_status, z_cvt, n, power_cvt, $
                              s_cvt, $
                              inverse_cvt, $
                              result
   ; Return
   RETURN, result
END
