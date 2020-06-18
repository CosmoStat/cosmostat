; $Id: //depot/idl/IDL_70/idl_src/libs/imsl/imsl_6.0/lib/imsl_covariances.pro#1 $      
;
; Copyright (c) 1970-2006, VISUAL NUMERICS Inc. All Rights Reserved.
; This software is confidential information which is proprietary to and a
; trade secret of Visual Numerics Inc.  Use, duplication or disclosure is
; subject to the terms of an appropriate license agreement.
;
FUNCTION imsl_covariances, X, $             ;INPUT 2-D array: floating point 
           corrected_sscp=corrected_sscp, $ ;INPUT Scalar ON/OFF flag
           correlation=correlation, $       ;INPUT Scalar ON/OFF flag
           frequencies=frequencies, $       ;INPUT 1-D array: floating point 
           incidence_mat=incidence_mat, $   ;OUTPUT 2-D array: floating point 
           means=means, $                   ;OUTPUT 1-D array: floating point 
           missing_val=missing_val, $       ;INPUT Scalar LONG
           sum_weights=sum_weights, $       ;OUTPUT Scalar floating point
           nmissing=nmissing, $             ;OUTPUT Scalar LONG
           nobs=nobs, $                     ;OUTPUT Scalar LONG
           stdev_correlation=stdev_correlation, $ ;INPUT Scalar ON/OFF flag
           var_covar=var_covar, $           ;INPUT Scalar ON/OFF flag
           weights=weights, $               ;INPUT 1-D array: floating point 
           double=double                    ;INPUT Scalar ON/OFF flag
           
@imsl_init.pro
   ON_ERROR, on_err_action
   ;
   ; Error checking:
   ;  - X must be either a 1-D or 2-D array. (m x n)
   ;  - If WEIGHTS is set, then it must be a 1-D array of length m.
   ;  - If FREQUENCIES is set, then it must be a 1-D array of length m.
   ;  - At most one of VAR_COVAR, CORRECTED_SSCP, CORRELATION, or
   ;    STD_CORRELATION can be specified.
   ;          
   nargs = n_params()
   IF (nargs NE 1)  THEN message, 'Incorrect number of arguments.'
     
   size_x = IMSL_LONG(size(x))
   IF ((size_x(0) NE 1) AND (size_x(0) NE 2)) THEN $
     message, 'X must be a 1-D or 2-D array.'
   m = IMSL_LONG(size_x(1))
   IF (size_x(0) EQ 2) THEN n = IMSL_LONG(size_x(2)) ELSE n = IMSL_1
   IF (KEYWORD_SET(weights)) THEN BEGIN 
      size_weights = IMSL_SIZE(weights)
      IF (size_weights(0) NE 1) THEN message, 'WEIGHTS must be a 1-D array.'
      IF (size_weights(1) NE m) THEN $
        message, 'WEIGHTS is not the correct size.'
   END
   IF (KEYWORD_SET(frequencies)) THEN BEGIN 
      size_freq = IMSL_LONG(size(frequencies))
      IF (size_freq(0) NE 1) THEN message, 'FREQUENCIES must be a 1-D array.'
      IF (size_freq(1) NE m) THEN $
        message, 'FREQUENCIES is not the correct size.'
   END
   IF ((KEYWORD_SET(VAR_COVAR) + KEYWORD_SET(CORRECTED_SSCP) + $
        KEYWORD_SET(CORRELATION) + KEYWORD_SET(STD_CORRELATION)) GT 1) THEN $
     message, 'At most one of VAR_COVAR, CORRECTED_SSCP, CORRELATION, OR '+ $
              'STD_CORRELATION can be specified.'
   ;
   ; Decide on what precision to use.
   ;
   type = TYP_FLOAT
   IF (size_x(N_ELEMENTS(size_x)-2) EQ  TYP_DOUBLE) THEN type = TYP_DOUBLE
   IF (KEYWORD_SET(double) EQ true) THEN type = TYP_DOUBLE
   ;
   ; Setup the parameters for the call to the system function.
   ;
   IF (ARG_PRESENT(nobs)) THEN nobs_spc = IMSL_0
   IF (ARG_PRESENT(nmissing)) THEN nmissing_spc = IMSL_0
   IF (KEYWORD_SET(missing_val)) THEN missing_val_cvt = IMSL_LONG(missing_val(0))
   IF (KEYWORD_SET(corrected_sscp)) THEN cor_sscp_cvt = IMSL_1
   IF (KEYWORD_SET(correlation)) THEN correlation_cvt = IMSL_1
   IF (KEYWORD_SET(stdev_correlation)) THEN stdev_corr_cvt = IMSL_1
   IF (KEYWORD_SET(var_covar)) THEN var_covar_cvt = IMSL_1
   ; Floating point arguments and keywords
   ;
   IF (type EQ TYP_DOUBLE) THEN BEGIN
      result = dblarr(n, n)
      x_cvt = double(transpose(x))
      IF (KEYWORD_SET(frequencies)) THEN frequencies_cvt = double(frequencies)
      IF (KEYWORD_SET(weights)) THEN weights_cvt = double(weights)
      IF (ARG_PRESENT(means)) THEN means_spc = dblarr(n)
      IF (ARG_PRESENT(sum_weights)) THEN sw_spc = DOUBLE(0)
      IF (ARG_PRESENT(incidence_mat)) THEN inc_mat_spc = dblarr(n, n)
   END ELSE BEGIN
      result = fltarr(n, n)
      x_cvt = float(transpose(x))
      IF (KEYWORD_SET(frequencies)) THEN frequencies_cvt = float(frequencies)
      IF (KEYWORD_SET(weights)) THEN weights_cvt = float(weights)
      IF (ARG_PRESENT(means)) THEN means_spc = fltarr(n)
      IF (ARG_PRESENT(sum_weights)) THEN sw_spc = FLOAT(0)
      IF (ARG_PRESENT(incidence_mat)) THEN inc_mat_spc = fltarr(n, n)
   END
   ;
   ; Call the system function.
   ;
   err_status = 0L
   MATHSTAT_126, type, err_status, X_cvt, m, n, $
                              correlation_cvt, $
                              frequencies_cvt, $
                              inc_mat_spc, $
                              means_spc, $
                              missing_val_cvt, $
                              nobs_spc, $
                              cor_sscp_cvt, $
                              stdev_corr_cvt, $
                              var_covar_cvt, $
                              weights_cvt, $
                              sw_spc, $
                              nmissing_spc, $
                              result

   IF (ARG_PRESENT(nobs)) THEN nobs = nobs_spc
   IF (ARG_PRESENT(sum_weights)) THEN sum_weights = sw_spc
   IF (ARG_PRESENT(nmissing)) THEN nmissing = nmissing_spc
   IF (ARG_PRESENT(means)) THEN means = means_spc
   IF (ARG_PRESENT(incidence_mat)) THEN incidence_mat = inc_mat_spc

; Return
   RETURN, result
END
