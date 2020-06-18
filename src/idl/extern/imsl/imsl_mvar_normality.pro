; $Id: //depot/idl/IDL_70/idl_src/libs/imsl/imsl_6.0/lib/imsl_mvar_normality.pro#1 $
;
; Copyright (c) 1970-2006, VISUAL NUMERICS Inc. All Rights Reserved.
; This software is confidential information which is proprietary to and a
; trade secret of Visual Numerics Inc.  Use, duplication or disclosure is
; subject to the terms of an appropriate license agreement.
;
FUNCTION imsl_mvar_normality, X, $          ;INPUT 2-D array: floating point 
           frequencies=frequencies, $       ;INPUT 1-D array: floating point 
           weights=weights, $               ;INPUT 1-D array: floating point 
           sum_weights=sum_weights, $       ;OUTPUT Scalar floating point
           sum_freqs=sum_freqs, $           ;OUTPUT Scalar floating point
           nmissing=nmissing, $             ;OUTPUT Scalar LONG
           means=means, $                   ;OUTPUT 1-D array: floating point 
           r_matrix=r_matrix, $             ;OUTPUT 2-D array: floating point 
           double=double                    ;INPUT Scalar ON/OFF flag
           
@imsl_init.pro
   ON_ERROR, on_err_action
   ;
   ; Error checking:
   ;  - X must be either a 1-D or 2-D array. (m x n)
   ;  - If WEIGHTS is set, then it must be a 1-D array of length m.
   ;  - If FREQUENCIES is set, then it must be a 1-D array of length m.
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
   ;
   ; Decide on what precision to use.
   ;
   type = TYP_FLOAT
   IF (size_x(N_ELEMENTS(size_x)-2) EQ  TYP_DOUBLE) THEN type = TYP_DOUBLE
   IF (KEYWORD_SET(double) EQ true) THEN type = TYP_DOUBLE
   ;
   ; Setup the parameters for the call to the system function.
   ;
   IF (ARG_PRESENT(nmissing)) THEN nmissing_spc = IMSL_0
   IF (ARG_PRESENT(sum_freqs)) THEN sf_spc = IMSL_0
   ; Floating point arguments and keywords
   ;
   IF (type EQ TYP_DOUBLE) THEN BEGIN
      result = dblarr(13)
      x_cvt = double(transpose(x))
      IF (KEYWORD_SET(weights)) THEN weights_cvt = double(weights)
      IF (KEYWORD_SET(frequencies)) THEN frequencies_cvt = double(frequencies)
      IF (ARG_PRESENT(means)) THEN means_spc = dblarr(n)
      IF (ARG_PRESENT(sum_weights)) THEN sw_spc = DOUBLE(0)
      IF (ARG_PRESENT(r_matrix)) THEN r_mat_spc = dblarr(n, n)
   END ELSE BEGIN
      result = fltarr(13)
      x_cvt = float(transpose(x))
      IF (KEYWORD_SET(weights)) THEN weights_cvt = float(weights)
      IF (KEYWORD_SET(frequencies)) THEN frequencies_cvt = float(frequencies)
      IF (ARG_PRESENT(means)) THEN means_spc = fltarr(n)
      IF (ARG_PRESENT(sum_weights)) THEN sw_spc = float(0)
      IF (ARG_PRESENT(r_matrix)) THEN r_mat_spc = fltarr(n, n)
   END
   ;
   ; Call the system function.
   ;
   err_status = 0L
   MATHSTAT_270, type, err_status, X_cvt, m, n, $
                              weights_cvt, $
                              frequencies_cvt, $
                              r_mat_spc, $
                              means_spc, $
                              sw_spc, $
                              sf_spc, $
                              nmissing_spc, $
                              result

   IF (ARG_PRESENT(sum_weights)) THEN sum_weights = sw_spc
   IF (ARG_PRESENT(sum_freqs)) THEN sum_freqs = sf_spc
   IF (ARG_PRESENT(nmissing)) THEN nmissing = nmissing_spc
   IF (ARG_PRESENT(means)) THEN means  =  means_spc   
   IF (ARG_PRESENT(r_matrix)) THEN r_matrix = TRANSPOSE(r_mat_spc)

; Return
   RETURN, result
END
