; $Id: //depot/idl/IDL_70/idl_src/libs/imsl/imsl_6.0/lib/imsl_chisqtest.pro#1 $
;
; Copyright (c) 1970-2006, VISUAL NUMERICS Inc. All Rights Reserved.
; This software is confidential information which is proprietary to and a
; trade secret of Visual Numerics Inc.  Use, duplication or disclosure is
; subject to the terms of an appropriate license agreement.
;
FUNCTION imsl_chisqtest, f, $                 ;INPUT Scalar STRING
                n_categories, $            ;INPUT Scalar LONG
                X, $                       ;INPUT 1-D array: floating point
                double=double, $           ;INPUT Scalar ON/OFF flag
                cutpoints=cutpoints, $     ;INPUT 1-D array: floating point
                frequencies=frequencies, $ ;INPUT 1-D array: floating point
                lower_bound=lower_bound, $ ;INPUT Scalar floating point
                n_params_estimated=n_params_estimated, $ ;INPUT Scalar LONG
                upper_bound=upper_bound, $ ;INPUT Scalar floating point
                equal_cutpoints=equal_cutpoints, $ ;INPUT Scalar ON/OFF flag
                cell_chisq=cell_chisq, $   ;OUTPUT 1-D array: floating point
                cell_counts=cell_counts, $ ;OUTPUT 1-D array: floating point
                cell_expected=cell_expected, $ ;OUTPUT 1-D array: floating point
                chi_squared=chi_squared, $ ;OUTPUT Scalar floating point
                df=df, $                   ;OUTPUT Scalar floating point
                used_cutpoints=used_cutpoints ;OUTPUT 1-D array: floating point

@imsl_init.pro
   ON_ERROR, on_err_action
   ;
   ; Error checking.
   ; - F must be a scalar string.
   ; - n_categories is converted to a scalar long. 
   ; - X must be a 1-D array. (set nobs = length)
   ; - Keywords EQUAL_CUTPOINTS and CUTPOINTS are mutually exclusive.
   ; - if CUTPOINTS is supplied, it must be a 1-D array of 
   ;   length (n_categories-1)
   ; - if FREQUENCIES is supplied, it must be a 1-D array of 
   ;   length (nobs)
   ; - The keywords UPPERBOUND and LOWERBOUND must be supplied together.
   ;                       
   nargs = n_params()
   IF (nargs NE 3) THEN message, 'Incorrect number of arguments.'
   size_f = IMSL_SIZE(f)
   IF ((N_ELEMENTS(f) NE 1) OR (size_f(N_ELEMENTS(size_f)-2) NE 7)) THEN $
     message, 'F must be a scalar string.'
   n_cat_cvt = IMSL_LONG(n_categories(0))
   size_x = IMSL_LONG(size(x))
   IF (size_x(0) NE 1) THEN message, 'X must be a 1-D array.'
   nobs = IMSL_LONG(size_x(1))
   IF ((KEYWORD_SET(cutpoints) + KEYWORD_SET(equal_cutpoints)) EQ 2) THEN $
     message, 'CUTPOINTS and EQUAL_CUTPOINTS are mutually exclusive.'
   IF (KEYWORD_SET(cutpoints)) THEN BEGIN 
      size_ctpts = IMSL_LONG(size(cutpoints))
      IF (size_ctpts(0) NE 1) THEN message, 'CUTPOINTS must be a 1-D array.'
      IF (size_ctpts(1) NE (n_cat_cvt-1)) THEN $
        message, 'CUTPOINTS is not the correct size.'
   END
   IF (KEYWORD_SET(frequencies)) THEN BEGIN 
      size_freq = IMSL_LONG(size(frequencies))
      IF (size_freq(0) NE 1) THEN message, 'FREQUENCIES must be a 1-D array.'
      IF (size_freq(1) NE nobs) THEN $
        message, 'FREQUENCIES is not the correct size.'
   END
   IF ((KEYWORD_SET(upper_bound) + KEYWORD_SET(lower_bound)) EQ 1) THEN $
     message, 'LOWER_BOUND AND UPPER_BOUND must be used together.'
   ; 
   ; Decide on what precision to use.
   type = TYP_FLOAT
   IF (size_x(N_ELEMENTS(size_x)-2) EQ  TYP_DOUBLE) THEN type = TYP_DOUBLE
   IF (KEYWORD_SET(double) EQ true) THEN type = TYP_DOUBLE
   ;
   ; Setup the parameters for the call to the system function.
   ; Input LONG keyword(s)
   ;
                used_cutpoints=used_cutpoints ;OUTPUT 1-D array: floating point
   IF (KEYWORD_SET(n_params_estimated)) THEN n_params_cvt = IMSL_LONG(n_params_estimated(0))
   IF (KEYWORD_SET(equal_cutpoints)) THEN eq_ctpts_cvt = IMSL_1
   ; Floating point arguments and keywords   
   IF (type EQ TYP_DOUBLE) THEN BEGIN
      result = double(0.0)
      x_cvt = double(x)
      IF (KEYWORD_SET(lower_bound)) THEN lbound_cvt = double(lower_bound)
      IF (KEYWORD_SET(upper_bound)) THEN ubound_cvt = double(upper_bound)
      IF (KEYWORD_SET(frequencies)) THEN frequencies_cvt = double(frequencies)
      IF (ARG_PRESENT(chi_squared)) THEN chi_squared_spc = double(0.0)
      IF (ARG_PRESENT(df)) THEN df_spc = double(0.0)
      IF (ARG_PRESENT(cell_chisq)) THEN cell_chisq_spc = dblarr(n_cat_cvt)
      IF (ARG_PRESENT(cell_counts)) THEN cell_counts_spc = dblarr(n_cat_cvt)
      IF (ARG_PRESENT(cell_expected)) THEN cell_expected_spc = dblarr(n_cat_cvt)
      ;Always send used_cutpoints.
      used_cutpoints_spc = dblarr(n_cat_cvt-1)
      IF (KEYWORD_SET(cutpoints)) THEN used_cutpoints_spc(*) = double(cutpoints) $
        ELSE eq_ctpts_cvt = IMSL_1
   END ELSE BEGIN
      result = float(0.0)
      x_cvt = float(x)
      IF (KEYWORD_SET(lower_bound)) THEN lbound_cvt = float(lower_bound)
      IF (KEYWORD_SET(upper_bound)) THEN ubound_cvt = float(upper_bound)
      IF (KEYWORD_SET(frequencies)) THEN frequencies_cvt = float(frequencies)
      IF (ARG_PRESENT(chi_squared)) THEN chi_squared_spc = float(0.0)
      IF (ARG_PRESENT(df)) THEN df_spc = float(0.0)
      IF (ARG_PRESENT(cell_chisq)) THEN cell_chisq_spc = fltarr(n_cat_cvt)
      IF (ARG_PRESENT(cell_counts)) THEN cell_counts_spc = fltarr(n_cat_cvt)
      IF (ARG_PRESENT(cell_expected)) THEN cell_expected_spc = fltarr(n_cat_cvt)
      ;Always send used_cutpoints.
      used_cutpoints_spc = fltarr(n_cat_cvt-1)
      IF (KEYWORD_SET(cutpoints)) THEN used_cutpoints_spc(*) = float(cutpoints) $
        ELSE eq_ctpts_cvt = IMSL_1
   END
   ;
   ; Call the system function.
   ;
   err_status = 0L
   MATHSTAT_118, type, err_status, $
                              f, $
                              n_cat_cvt, $
                              X_cvt, $
                              nobs, $
                              frequencies_cvt, $
                              lbound_cvt, $
                              n_params_cvt, $
                              ubound_cvt, $
                              eq_ctpts_cvt, $
                              cell_chisq_spc, $
                              cell_counts_spc, $
                              cell_expected_spc, $
                              chi_squared_spc, $
                              df_spc, $
                              used_cutpoints_spc, $
                              result
   
   IF (ARG_PRESENT(chi_squared)) THEN chi_squared = chi_squared_spc
   IF (ARG_PRESENT(df)) THEN df = df_spc
   IF (ARG_PRESENT(cell_chisq)) THEN cell_chisq = cell_chisq_spc
   IF (ARG_PRESENT(cell_counts)) THEN cell_counts = cell_counts_spc
   IF (ARG_PRESENT(cell_expected)) THEN cell_expected = cell_expected_spc
   IF (ARG_PRESENT(used_cutpoints)) THEN used_cutpoints = used_cutpoints_spc

   RETURN, result
END

