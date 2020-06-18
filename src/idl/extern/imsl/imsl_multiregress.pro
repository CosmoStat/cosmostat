; $Id: //depot/idl/IDL_70/idl_src/libs/imsl/imsl_6.0/lib/imsl_multiregress.pro#1 $      
;
; Copyright (c) 1970-2006, VISUAL NUMERICS Inc. All Rights Reserved.
; This software is confidential information which is proprietary to and a
; trade secret of Visual Numerics Inc.  Use, duplication or disclosure is
; subject to the terms of an appropriate license agreement.
;
FUNCTION imsl_multiregress, X, $                ;INPUT 1-D array: floating point 
             Y, $                           ;INPUT 1-D or 2-D array: floating point 
             weights=weights, $             ;INPUT 1-D array: floating point 
             tolerance=tolerance, $         ;INPUT Scalar floating point 
             frequencies=frequencies, $     ;INPUT 1-D array: floating point 
             no_intercept=no_intercept, $   ;INPUT Scalar ON/OFF flag
             residual=residual, $           ;OUTPUT 1-D array: floating point 
             anova_table=anova_table, $     ;OUTPUT 1-D array: floating point 
             predict_info=predict_info, $   ;OUTPUT 1-D array: BYTE
             xmean=xmean, $                 ;OUTPUT Scalar floating point 
             coef_covariances=coef_covariances, $ ;OUTPUT 2-D array: floating point 
             coef_vif=coef_vif, $           ;OUTPUT 1-D array: floating point 
             rank=rank, $                   ;OUTPUT Scalar LONG
             t_tests=t_tests, $             ;OUTPUT 1-D array: floating point 
             double=double                  ;INPUT Scalar ON/OFF flag

@imsl_init.pro
   ON_ERROR, on_err_action
   ;
   ; Error checking:
   ;  - X must be either a 1-D or 2-D array. (m x n)
   ;  - Y must be either a 1-D or 2-D array. (m x ndep)
   ;  - If WEIGHTS is set, then it must be a 1-D array of length m.
   ;  - If FREQUENCIES is set, then it must be a 1-D array of length m.
   ;          
   nargs = n_params()
   IF (nargs NE 2)  THEN message, 'Incorrect number of arguments.'
     
   size_x = IMSL_LONG(size(x))
   IF ((size_x(0) NE 1) AND (size_x(0) NE 2)) THEN $
     message, 'X must be a 1-D or 2-D array.'
   m = IMSL_LONG(size_x(1))
   IF (size_x(0) EQ 2) THEN n = IMSL_LONG(size_x(2)) ELSE n = IMSL_1
   size_y = IMSL_LONG(size(y))
   IF ((size_y(0) NE 1) AND (size_y(0) NE 2)) THEN $
     MESSAGE,  'Y must be a 1-D or 2-D array.'
   IF (size_y(1) NE m) THEN MESSAGE,  'Y is not the correct size.'
   IF (SIZE_y(0) EQ 2) THEN ndep = size_y(2) ELSE ndep = IMSL_1
   IF (KEYWORD_SET(weights)) THEN BEGIN 
      size_weights = IMSL_SIZE(weights)
      IF (size_weights(0) NE 1) THEN message, 'WEIGHTS must be a 1-D array.'
      IF (size_weights(1) NE m) THEN $
        message, 'WEIGHTS is not the correct size.'
   END
   IF (KEYWORD_SET(frequencies)) THEN BEGIN 
      size_frequencies = IMSL_SIZE(frequencies)
      IF (size_frequencies(0) NE 1) THEN message, 'FREQUENCIES must be a 1-D array.'
      IF (size_frequencies(1) NE m) THEN $
        message, 'FREQUENCIES is not the correct size.'
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
   IF (ARG_PRESENT(rank)) THEN rank_spc = IMSL_0
   IF (ARG_PRESENT(predict_info)) THEN predict_spc = byte(0)
   intcpt_term = IMSL_1
   IF (KEYWORD_SET(no_intercept)) THEN BEGIN
      no_intrcpt_cvt = IMSL_1
      intcpt_term = IMSL_0
   END
   
   ; Floating point arguments and keywords
   ;
   IF (type EQ TYP_DOUBLE) THEN BEGIN
      result = dblarr(n + intcpt_term, ndep)
      x_cvt = double(transpose(x))
      IF (ndep EQ 2) THEN y_cvt = double(TRANSPOSE(y)) ELSE y_cvt = double(y)
      IF (KEYWORD_SET(weights)) THEN weights_cvt = double(weights)
      IF (KEYWORD_SET(frequencies)) THEN frequencies_cvt = double(frequencies)
      IF (KEYWORD_SET(tolerance)) THEN tolerance_cvt = double(tolerance(0))
      IF (ARG_PRESENT(anova_table)) THEN anova_table_spc = dblarr(ndep, 15)
      IF (ARG_PRESENT(residual)) THEN residual_spc = dblarr(m)
      IF (ARG_PRESENT(coef_covariances)) THEN $
        coef_cov_spc = dblarr(n+intcpt_term, n+intcpt_term, ndep)
      IF (ARG_PRESENT(coef_vif)) THEN coef_vif_spc = dblarr(n+intcpt_term)
      IF (ARG_PRESENT(t_tests)) THEN t_tests_spc = dblarr(4, n+intcpt_term)
   END ELSE BEGIN
      result = fltarr(n + intcpt_term, ndep)
      x_cvt = float(transpose(x))
      IF (ndep EQ 2) THEN y_cvt = float(TRANSPOSE(y)) else y_cvt = float(y)
      IF (KEYWORD_SET(weights)) THEN weights_cvt = float(weights)
      IF (KEYWORD_SET(frequencies)) THEN frequencies_cvt = float(frequencies)
      IF (KEYWORD_SET(tolerance)) THEN tolerance_cvt = float(tolerance(0))
      IF (ARG_PRESENT(anova_table)) THEN anova_table_spc = fltarr(ndep, 15)
      IF (ARG_PRESENT(residual)) THEN residual_spc = fltarr(m)
      IF (ARG_PRESENT(coef_covariances)) THEN $
        coef_cov_spc = fltarr(n+intcpt_term, n+intcpt_term, ndep)
      IF (ARG_PRESENT(coef_vif)) THEN coef_vif_spc = fltarr(n+intcpt_term)
      IF (ARG_PRESENT(t_tests)) THEN t_tests_spc = fltarr(4, n+intcpt_term)
   END
   ;
   ; Call the system function.
   ;
   err_status = 0L
   MATHSTAT_162, type, err_status, X_cvt, $
                              Y_cvt, $
                              m, $
                              n, $
                              ndep, $
                              weights_cvt, $
                              tolerance_cvt, $
                              frequencies_cvt, $
                              no_intrcpt_cvt, $
                              residual_spc, $
                              anova_table_spc, $
                              predict_spc, $
                              xmean_spc, $
                              coef_cov_spc, $
                              coef_vif_spc, $
                              rank_spc, $
                              t_tests_spc, $
                              result

   IF (ARG_PRESENT(predict_info)) THEN predict_info = predict_spc
   IF (ARG_PRESENT(anova_table)) THEN anova_table = TRANSPOSE(anova_table_spc)
   IF (ARG_PRESENT(residual)) THEN residual = residual_spc
   IF (ARG_PRESENT(xmean)) THEN xmean = xmean_spc
   IF (ARG_PRESENT(rank)) THEN rank = rank_spc
   IF (ARG_PRESENT(coef_vif)) THEN coef_vif = coef_vif_spc   
   IF (ARG_PRESENT(coef_covariances)) THEN coef_covariances = coef_cov_spc   
   IF (ARG_PRESENT(t_tests)) THEN t_tests = transpose(t_tests_spc)

   IF (ndep GT 1) THEN result = TRANSPOSE(result)
; Return
   RETURN, result
END
