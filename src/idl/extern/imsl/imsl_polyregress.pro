; $Id: //depot/idl/IDL_70/idl_src/libs/imsl/imsl_6.0/lib/imsl_polyregress.pro#1 $      
;
; Copyright (c) 1970-2006, VISUAL NUMERICS Inc. All Rights Reserved.
; This software is confidential information which is proprietary to and a
; trade secret of Visual Numerics Inc.  Use, duplication or disclosure is
; subject to the terms of an appropriate license agreement.
;
FUNCTION imsl_polyregress, X, $                ;INPUT 1-D array: floating point 
             Y, $                           ;INPUT 1-D array: floating point 
             Degree, $                      ;INPUT Scalar LONG
             weights=weights, $             ;INPUT 1-D array: 
             residual=residual, $           ;OUTPUT 1-D array: floating point 
             anova_table=anova_table, $     ;OUTPUT 1-D array: floating point 
             df_pure_error=df_pure_error, $ ;OUTPUT Scalar LONG
             predict_info=predict_info, $   ;OUTPUT 1-D array: BYTE
             ssq_lof=ssq_lof, $             ;OUTPUT 2-D array: floating point 
             ssq_poly=ssq_poly, $           ;OUTPUT 2-D array: floating point 
             ssq_pure_error=ssq_pure_error,$ ;OUTPUT Scalar floating point 
             xmean=xmean, $                 ;OUTPUT Scalar floating point 
             xvariance=xvariance, $         ;OUTPUT Scalar floating point 
             double=double                  ;INPUT Scalar ON/OFF flag
@imsl_init.pro
   ON_ERROR, on_err_action
   ;
   ; Error checking:
   ;  - X must be a 1-D array. (set n = length)
   ;  - Y must be a 1-D array of length = n.
   ;  - Degree must be a scalar greater than or equal to 0.
   ;  - If WEIGHTS is set, then it must be a 1-D array of length m.
   ;          
   nargs = n_params()
   IF (nargs NE 3)  THEN message, 'Incorrect number of arguments.'
     
   size_x = IMSL_LONG(size(x))
   IF (size_x(0) NE 1) THEN $
     message, 'X must be a 1-D array.'
   n = IMSL_LONG(size_x(1))
   size_y = IMSL_LONG(size(y))
   IF (size_y(0) NE 1) THEN $
     message, 'Y must be a 1-D array.'
   IF (N_ELEMENTS(y) NE n) THEN message, 'Y is not the correct size.'
   degree_cvt = IMSL_LONG(degree(0))
   IF (degree_cvt LT 0) THEN message,'DEGREE must be positive.'
   IF (KEYWORD_SET(weights)) THEN BEGIN 
      size_weights = IMSL_SIZE(weights)
      IF (size_weights(0) NE 1) THEN message, 'WEIGHTS must be a 1-D array.'
      IF (size_weights(1) NE n) THEN $
        message, 'WEIGHTS is not the correct size.'
   END
   ;
   ; Decide on what precision to use.
   ;
   type = TYP_FLOAT
   IF (size_x(N_ELEMENTS(size_x)-2) EQ  TYP_DOUBLE) THEN type = TYP_DOUBLE
   IF (size_y(N_ELEMENTS(size_y)-2) EQ  TYP_DOUBLE) THEN type = TYP_DOUBLE
   IF (KEYWORD_SET(double) EQ true) THEN type = TYP_DOUBLE
   ;
   ; Setup the parameters for the call to the system function.
   ;
   IF (ARG_PRESENT(df_pure_error)) THEN df_error_spc = IMSL_0
   IF (ARG_PRESENT(predict_info)) THEN predict_spc = byte(0)
   ; Floating point arguments and keywords
   ;
   IF (type EQ TYP_DOUBLE) THEN BEGIN
      result = dblarr(degree_cvt+1)
      x_cvt = double(x)
      y_cvt = double(y)
      IF (KEYWORD_SET(weights)) THEN weights_cvt = double(weights)
      IF (ARG_PRESENT(ssq_poly)) THEN ssq_poly_spc = dblarr(4, degree_cvt)
      IF (ARG_PRESENT(ssq_lof)) THEN ssq_lof_spc = dblarr(4, degree_cvt)
      IF (ARG_PRESENT(anova_table)) THEN anova_table_spc = dblarr(15)
      IF (ARG_PRESENT(residual)) THEN residual_spc = dblarr(n)
      IF (ARG_PRESENT(xmean)) THEN xmean_spc = double(0.0)
      IF (ARG_PRESENT(xvariance)) THEN xvariance_spc = double(0.0)
      IF (ARG_PRESENT(ssq_pure_error)) THEN ssq_error_spc = double(0.0)
   END ELSE BEGIN
      result = fltarr(degree_cvt+1)
      x_cvt = float(x)
      y_cvt = float(y)
      IF (KEYWORD_SET(weights)) THEN weights_cvt = float(weights)
      IF (ARG_PRESENT(ssq_poly)) THEN ssq_poly_spc = fltarr(4, degree_cvt)
      IF (ARG_PRESENT(ssq_lof)) THEN ssq_lof_spc = fltarr(4, degree_cvt)
      IF (ARG_PRESENT(anova_table)) THEN anova_table_spc = fltarr(15)
      IF (ARG_PRESENT(residual)) THEN residual_spc = fltarr(n)
      IF (ARG_PRESENT(xmean)) THEN xmean_spc = float(0.0)
      IF (ARG_PRESENT(xvariance)) THEN xvariance_spc = float(0.0)
      IF (ARG_PRESENT(ssq_pure_error)) THEN ssq_error_spc = float(0.0)
   END
   ;
   ; Call the system function.
   ;
   err_status = 0L
   MATHSTAT_174, type, err_status, X_cvt, $
                              Y_cvt, $
                              Degree_cvt, $
                              n, $
                              weights_cvt, $
                              residual_spc, $
                              anova_table_spc, $
                              df_error_spc, $
                              predict_spc, $
                              ssq_lof_spc, $
                              ssq_poly_spc, $
                              ssq_error_spc,$
                              xmean_spc, $
                              xvariance_spc, $
                              result

   IF (ARG_PRESENT(df_pure_error)) THEN df_pure_error = df_error_spc
   IF (ARG_PRESENT(predict_info)) THEN predict_info = predict_spc
   IF (ARG_PRESENT(anova_table)) THEN anova_table = anova_table_spc
   IF (ARG_PRESENT(residual)) THEN residual = residual_spc
   IF (ARG_PRESENT(xmean)) THEN xmean = xmean_spc
   IF (ARG_PRESENT(xvariance)) THEN xvariance = xvariance_spc
   IF (ARG_PRESENT(ssq_pure_error)) THEN ssq_pure_error = ssq_error_spc
   IF (ARG_PRESENT(ssq_poly)) THEN ssq_poly = transpose(ssq_poly_spc)
   IF (ARG_PRESENT(ssq_lof)) THEN ssq_lof = transpose(ssq_lof_spc)
   
   
; Return
   RETURN, result
END
