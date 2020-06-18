; $Id: //depot/idl/IDL_70/idl_src/libs/imsl/imsl_6.0/lib/imsl_simplestat.pro#1 $      
;
; Copyright (c) 1970-2006, VISUAL NUMERICS Inc. All Rights Reserved.
; This software is confidential information which is proprietary to and a
; trade secret of Visual Numerics Inc.  Use, duplication or disclosure is
; subject to the terms of an appropriate license agreement.
;
FUNCTION imsl_simplestat, X, $                    ;INPUT 1-D or 2-D array: floating point 
           double=double, $                    ;INPUT Scalar ON/OFF flag
           conf_means=conf_means, $            ;INPUT Scalar floating point
           conf_variances=conf_variances, $    ;INPUT Scalar floating point
           elementwise=elementwise, $          ;INPUT Scalar ON/OFF flag
           frequencies=frequencies, $          ;INPUT 1-D array: floating point
           median_and_scale=median_and_scale, $ ;INPUT Scalar ON/OFF flag
           median_only=median_only, $          ;INPUT Scalar ON/OFF flag
           weights=weights                     ;INPUT 1-D array: floating point

@imsl_init.pro
   ON_ERROR, on_err_action
   ;
   ; Error checking:
   ;  - X must be either a 1-D or 2-D array. (m x n)
   ;  - MEDIAN_ONLY and MEDIAN_AND_SCALE are mutually exclusive.
   ;  - If WEIGHTS is set, then it must be a 1-D array of length m.
   ;  - If FREQUENCIES is set, then it must be a 1-D array of length m.
   ;          
   nargs = n_params()
   IF (nargs NE 1)  THEN message, 'Incorrect number of arguments.'
     
   size_x = IMSL_LONG(size(x))
   ndim = IMSL_LONG(size_x(0))
   IF ((ndim NE 1) AND (ndim NE 2)) THEN message, 'X must be a 1-D or 2-D array.'
   m = IMSL_LONG(size_x(1))
   IF (ndim EQ 2) THEN n = IMSL_LONG(size_x(2)) ELSE n = IMSL_1
   IF ((KEYWORD_SET(median_only) + KEYWORD_SET(median_and_scale)) GT 1) THEN $
     message, 'Keywords MEDIAN_ONLY AND MEDIAN_AND_SCALE are mutually exclusive.'
   IF (KEYWORD_SET(weights)) THEN BEGIN 
      size_weights = IMSL_SIZE(weights)
      IF ((size_weights(0) NE 1) OR (N_ELEMENTS(weights) NE m)) THEN $
        message, 'WEIGHTS is not the correct size.'
   END
   IF (KEYWORD_SET(frequencies)) THEN BEGIN 
      size_frequencies = IMSL_SIZE(frequencies)
      IF ((size_frequencies(0) NE 1) OR (N_ELEMENTS(frequencies) NE m)) THEN $
        message, 'FREQUENCIES is not the correct size.'
   END
   ;
   ; Decide on what precision to use.
   ;
   type = TYP_FLOAT
   IF (size_x(N_ELEMENTS(size_x)-2) EQ  TYP_DOUBLE) THEN type = TYP_DOUBLE
   IF (KEYWORD_SET(double) EQ true) THEN type = TYP_DOUBLE
   ;
   ; Setup the parameters for the call to the system function.
   ;
   ; Input LONG
   IF (KEYWORD_SET(elementwise)) THEN elementwise_cvt = IMSL_1
   IF (KEYWORD_SET(median_and_scale)) THEN median_scl_cvt = IMSL_1
   IF (KEYWORD_SET(median_only)) THEN median_only_cvt = IMSL_1
   ;
   ; Floating point arguments and keywords
   ;
   dim1 = 14
   IF (KEYWORD_SET(median_and_scale)) THEN dim1 = 17
   IF (KEYWORD_SET(median_only)) THEN dim1 = 15

   IF (type EQ TYP_DOUBLE) THEN BEGIN
      result = dblarr(n, dim1)
      x_cvt = double(transpose(x))
      IF (KEYWORD_SET(conf_means)) THEN conf_means_cvt = double(conf_means) $
        ELSE conf_means_cvt = double(95.)
      IF (KEYWORD_SET(conf_variances)) THEN conf_var_cvt = double(conf_variances) $
        ELSE conf_var_cvt = double(95.)
      IF (KEYWORD_SET(weights)) THEN weights_cvt = double(weights)
      IF (KEYWORD_SET(frequencies)) THEN frequencies_cvt = double(frequencies)
   END ELSE BEGIN
      result = fltarr(n, dim1)
      x_cvt = float(transpose(x))
      IF (KEYWORD_SET(conf_means)) THEN conf_means_cvt = float(conf_means) $
        ELSE conf_means_cvt = float(95.)
      IF (KEYWORD_SET(conf_variances)) THEN conf_var_cvt = float(conf_variances) $
        ELSE conf_var_cvt = float(95.)
      IF (KEYWORD_SET(weights)) THEN weights_cvt = float(weights)
      IF (KEYWORD_SET(frequencies)) THEN frequencies_cvt = float(frequencies)
   END
   ;
   ; Call the system function.
   ;
   err_status = 0L
   MATHSTAT_187, type, err_status, X_cvt, m, n, $
                              conf_means_cvt, $
                              conf_var_cvt, $
                              elementwise_cvt, $
                              frequencies_cvt, $
                              median_scl_cvt, $
                              median_only_cvt, $
                              weights_cvt, $
                              result
   ; Return
   RETURN, transpose(result)
END
