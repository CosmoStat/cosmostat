; $Id: //depot/idl/IDL_70/idl_src/libs/imsl/imsl_6.0/lib/imsl_tcdf.pro#1 $
;
; Copyright (c) 1970-2006, VISUAL NUMERICS Inc. All Rights Reserved.
; This software is confidential information which is proprietary to and a
; trade secret of Visual Numerics Inc.  Use, duplication or disclosure is
; subject to the terms of an appropriate license agreement.
;
FUNCTION imsl_tcdf, x, $         ;INPUT Scalar or 1-D array: floating point
                 y, $         ;INPUT Scalar or 1-D array: floating point
                 z, $         ;INPUT Scalar or 1-D array: floating point
                double=double, $  ;INPUT Scalar ON/OFF flag
                inverse=inverse   ;INPUT Scalar ON/OFF flag
@imsl_init.pro
   ON_ERROR,  on_err_action
   
   nargs = n_params()
   IF (nargs EQ 2) THEN begin
    result = imsl_call_spcl_fcn2(x, y, $
                            fcn_name = "TCDF", $
                            double = double, $
                            inverse = inverse)
   END ELSE BEGIN
      ; 
      ; Non-Central T CDF
      ; Note that this case will use scalars only.
      ;
      IF (nargs EQ 3) THEN BEGIN
         IF (N_ELEMENTS(x) EQ FALSE) THEN message, "T is undefined"
         IF (N_ELEMENTS(y) EQ FALSE) THEN message, "DF is undefined"
         IF (N_ELEMENTS(z) EQ FALSE) THEN message, "DELTA is undefined"
      END ELSE message, "Incorrect number of arguments."
      ;
      ; Setup the parameters for the call to the system function.
      ; 
      ; Input LONG argument(s)
      ;
      df_cvt   =   IMSL_LONG(y(0))
      IF (KEYWORD_SET(inverse)) THEN inv_set = IMSL_1
   
      IF (KEYWORD_SET(double)) THEN begin
         type = TYP_DOUBLE
         t_cvt = DOUBLE(x(0))
         delta_cvt = DOUBLE(z(0))
         result = DOUBLE(0.0)
      END ELSE BEGIN
         type = TYP_FLOAT
         t_cvt = FLOAT(x(0))
         delta_cvt = FLOAT(z(0))
         result = FLOAT(0.0)
      END
      ;
      ; Call the system function.
      ;
      err_status = 0L
      MATHSTAT_263,  type,  err_status,  t_cvt,  df_cvt,  $
            delta_cvt, inv_set, result
   END    
   ;
   ; Return.
   ;
   RETURN, result
 END
