; $Id: //depot/idl/IDL_70/idl_src/libs/imsl/imsl_6.0/lib/imsl_stepwise.pro#1 $
;
; Copyright (c) 1970-2006, VISUAL NUMERICS Inc. All Rights Reserved.
; This software is confidential information which is proprietary to and a
; trade secret of Visual Numerics Inc.  Use, duplication or disclosure is
; subject to the terms of an appropriate license agreement.
;
PRO imsl_stepwise, x, $                          ;INPUT 2-D array: floating point
                   y, $                          ;INPUT 1-D array: floating point
                   all_steps=all_steps, $        ;INPUT Scalar ON/OFF flag
                   first_step=first_step, $      ;INPUT Scalar ON/OFF flag
                   inter_step=inter_step, $      ;INPUT Scalar ON/OFF flag
                   last_step=last_step, $        ;INPUT Scalar ON/OFF flag
                   forward=forward, $            ;INPUT Scalar ON/OFF flag
                   backward=backward, $          ;INPUT Scalar ON/OFF flag
                   stepwise=stepwise, $          ;INPUT Scalar ON/OFF flag
                   double=double, $              ;INPUT Scalar ON/OFF flag
                   n_steps=n_steps, $            ;INPUT Scalar LONG
                   force=force, $                ;INPUT Scalar LONG
                   cov_nobs=cov_nobs, $          ;INPUT Scalar LONG
                   level=level, $                ;INPUT 1-D array: LONG
                   p_in=p_in, $                  ;INPUT Scalar floating point
                   p_out=p_out, $                ;INPUT Scalar floating point
                   tolerance=tolerance, $        ;INPUT Scalar floating point
                   frequencies=frequencies, $    ;INPUT 1-D array: floating point
                   weights=weights, $            ;INPUT 1-D array: floating point
                   cov_input=cov_input, $        ;INPUT 2-D array: floating point
                   iend=iend, $                  ;OUTPUT Scalar LONG
                   anova_table = anova_table, $  ;OUTPUT 1-D array: floating point
                   coef_t_tests=coef_t_tests, $  ;OUTPUT 2-D array: floating point
                   coef_vif=coef_vif, $          ;OUTPUT 1-D array: floating point
                   cov_swept=cov_swept, $        ;OUTPUT 1-D array: floating point
                   history=history, $            ;OUTPUT 1-D array: floating point
                   swept=swept                   ;OUTPUT 2-D array: floating point
                    
@imsl_init.pro
   ON_ERROR, on_err_action
   ;
   ; Error checking:
   ; The following checks are performed.
   ;    o   x:  must be either a 1-D array, in which case n_candidate
   ;            is one, or a 2-D array.
   ;    o   y:  must be a 1-D array of the same length as the first
   ;            dimension of p_argv[0].
   ;    o   The following keywords are mutually exclusive:
   ;        FIRST_STEP, INTRER_STEP, LAST_STEP, ALL_STEPS.
   ;    o   If WEIGHTS is present, then weights must be a one dimensional 
   ;        array of length n_rows.
   ;    o   If FREQUENCIES is present, then weights must be a one dimensional 
   ;        array of length n_rows.
   ;    o   The following keywords are mutually exclusive:
   ;        FORWARD, BACKWARD, STEPWISE.
   ;    o   If LEVEL is present, it must be a 1-D array of length
   ;        n_candidate + 1.
   ;    o   COV_INPUT and COV_NOBS must be used together.
   ;    o   If COV_INPUT is present, then it must be a 2-D square matrix
   ;        of order (n_candidate+1).
   ;
   nargs = n_params()
   IF (nargs NE 2) THEN $
         message, "Incorrect number of arguments."
   size_x = IMSL_SIZE(x)
   IF ((size_x(0) NE 1) and (size_x(0) NE 2)) THEN BEGIN
      message, "X must be a 1-D or 2-D array."
   END
   n_rows = size_x(1)
   IF (size_x(0) EQ 1) THEN n_candidate = 1 ELSE n_candidate = size_x(2)
   size_y = IMSL_SIZE(y)
   IF (size_y(0) NE 1) THEN BEGIN
      message, "Y must be a 1-D array."
   END
   IF (size_y(1) NE n_rows) THEN $
     message, "The input array, Y, is not the correct length"
   IF ((KEYWORD_SET(first_step)+KEYWORD_SET(inter_step)+ $
       KEYWORD_SET(last_step)+KEYWORD_SET(all_steps)) GT 1) THEN $
     message, "The following keywords are mutually exclusive:  " + $
               "FIRST_STEP, INTRER_STEP, LAST_STEP, AND ALL_STEPS"
   IF (KEYWORD_SET(weights) EQ TRUE) THEN BEGIN
     size_tmp = IMSL_SIZE(weights)
     IF (size_tmp(0) NE 1) THEN BEGIN
        message, "The weights array must be 1-D"
     END
     IF (size_tmp(1) NE n_rows) THEN BEGIN
        message, "The length of the weights array be equal to the first dimension of X"
     END
   END
   IF (KEYWORD_SET(frequencies) EQ TRUE) THEN BEGIN
     size_tmp = IMSL_SIZE(frequencies)
     IF (size_tmp(0) NE 1) THEN BEGIN
        message, "The frequencies array must be 1-D"
     END
     IF (size_tmp(1) NE n_rows) THEN BEGIN
        message, "The length of the frequencies array must be equal to the first dimension of X"
     END
   END
   IF ((KEYWORD_SET(forward) + KEYWORD_SET(backward)+ KEYWORD_SET(stepwise)) GT 1) THEN $
    message, "The following keywords are mutually exclusive: FORWARD, BACKWARD, and STEPWISE"
   IF (KEYWORD_SET(level) EQ TRUE) THEN BEGIN
     size_tmp = IMSL_SIZE(level)
     IF (size_tmp(0) NE 1) THEN BEGIN
        message, "LEVEL must be a 1-D array."
     END
     IF (size_tmp(1) NE n_candidate+1) THEN BEGIN
        message, "The length of the LEVEL array is incorrect"
     END
   END
   IF ((KEYWORD_SET(cov_input) + KEYWORD_SET(cov_nobs)) EQ 1) THEN $
     message, "COV_NOBS and COV_INPUT must be used together"
   IF (KEYWORD_SET(cov_input) EQ TRUE) THEN BEGIN
      size_tmp = IMSL_SIZE(cov_input)
      IF (size_tmp(0) NE 2) THEN $
        message, "The COV_INPUT array must be 2-D"
      IF ((size_tmp(1) NE (n_candidate+1)) OR ((size_tmp(2) NE (n_candidate+1)))) THEN $
          message,  "The COV_INPUT array is not the correct size"
   END
   ;
   ; Decide on what precision to use.
   ; Use the highest precision of the input argument(s).
   ;
   type = TYP_FLOAT
   IF (size_x(N_ELEMENTS(size_x)-2) GT TYP_FLOAT) THEN type = TYP_DOUBLE
   IF (size_y(N_ELEMENTS(size_y)-2) GT TYP_FLOAT) THEN type = TYP_DOUBLE
   IF (KEYWORD_SET(double) EQ true) THEN type = TYP_DOUBLE
   ;
   ; Setup the parameters for the call to the system function.
   ;
   ; Input LONG keyword(s)
   IF (KEYWORD_SET(all_steps) EQ TRUE) THEN $
     all_steps_cvt = IMSL_1
   IF (KEYWORD_SET(first_step) EQ TRUE) THEN $
     first_step_cvt = IMSL_1
   IF (KEYWORD_SET(inter_step) EQ TRUE) THEN $
     inter_step_cvt = IMSL_1
   IF (KEYWORD_SET(last_step) EQ TRUE) THEN $
     last_step_cvt = IMSL_1
   IF (KEYWORD_SET(forward) EQ TRUE) THEN $
     forward_cvt = IMSL_1
   IF (KEYWORD_SET(backward) EQ TRUE) THEN $
     backward_cvt = IMSL_1
   IF (KEYWORD_SET(stepwise) EQ TRUE) THEN $
     stepwise_cvt = IMSL_1
   IF (KEYWORD_SET(n_steps) EQ TRUE) THEN $
     n_steps_cvt = (IMSL_LONG(n_steps))(0)
   IF (KEYWORD_SET(force) EQ TRUE) THEN $
     force_cvt = (IMSL_LONG(force))(0)
   IF (KEYWORD_SET(cov_nobs) EQ TRUE) THEN $
     cov_nobs_cvt = (IMSL_LONG(cov_nobs))(0)
   IF (KEYWORD_SET(level) EQ TRUE) THEN $
     level_cvt = (IMSL_LONG(level))(0)
   ;
   ; Output LONG keyword(s)
   IF (ARG_PRESENT(iend) EQ TRUE) THEN $
     iend_spc = IMSL_0
   ;
   ; Floating point arguments and keywords
   IF (type EQ TYP_DOUBLE) THEN BEGIN
      ; 
      ; Input 
      x_cvt = double(transpose(x))
      y_cvt = double(y)
      IF (ARG_PRESENT(p_in) EQ TRUE) THEN $
        p_in_cvt = (double(p_in))(0)
      IF (ARG_PRESENT(p_out) EQ TRUE) THEN $
        p_out_cvt = (double(p_out))(0)
      IF (ARG_PRESENT(tolerance) EQ TRUE) THEN $
        tolerance_cvt = (double(tolerance))(0)
      IF (ARG_PRESENT(frequencies) EQ TRUE) THEN $
        frequencies_cvt = double(frequencies)
      IF (ARG_PRESENT(weights) EQ TRUE) THEN $
        weights_cvt = double(weights)
      IF (ARG_PRESENT(cov_input) EQ TRUE) THEN $
        cov_input_cvt = double(cov_input)
      ;
      ; Output keywords.
      ;
      IF (ARG_PRESENT(anova_table) EQ TRUE) THEN $
        anova_table_spc = dblarr(16)
      IF (ARG_PRESENT(coef_t_tests) EQ TRUE) THEN $
        coef_t_spc = dblarr(4, n_candidate)
      IF (ARG_PRESENT(coef_vif) EQ TRUE) THEN $
        coef_vif_spc = dblarr(n_candidate)
      IF (ARG_PRESENT(cov_swept) EQ TRUE) THEN $
        cov_swept_spc = dblarr(n_candidate+1, n_candidate+1)
      IF (ARG_PRESENT(history) EQ TRUE) THEN $
        history_spc = dblarr(n_candidate+1)
      IF (ARG_PRESENT(swept) EQ TRUE) THEN $
        swept_spc = dblarr(n_candidate+1)
   END ELSE BEGIN
      ; 
      ; Input 
      x_cvt = float(transpose(x))
      y_cvt = float(y)
      IF (KEYWORD_SET(p_in) EQ TRUE) THEN $
        p_in_cvt = (float(p_in))(0)
      IF (KEYWORD_SET(p_out) EQ TRUE) THEN $
        p_out_cvt = (float(p_out))(0)
      IF (KEYWORD_SET(tolerance) EQ TRUE) THEN $
        tolerance_cvt = (float(tolerance))(0)
      IF (KEYWORD_SET(frequencies) EQ TRUE) THEN $
        frequencies_cvt = float(frequencies)
      IF (KEYWORD_SET(weights) EQ TRUE) THEN $
        weights_cvt = float(weights)
      IF (KEYWORD_SET(cov_input) EQ TRUE) THEN $
        cov_input_cvt = float(cov_input)
      ;
      ; Output keywords.
      ;
      IF (ARG_PRESENT(anova_table) EQ TRUE) THEN $
        anova_table_spc = fltarr(16)
      IF (ARG_PRESENT(coef_t_tests) EQ TRUE) THEN $
        coef_t_spc = fltarr(4, n_candidate)
      IF (ARG_PRESENT(coef_vif) EQ TRUE) THEN $
        coef_vif_spc = fltarr(n_candidate)
      IF (ARG_PRESENT(cov_swept) EQ TRUE) THEN $
        cov_swept_spc = fltarr(n_candidate+1, n_candidate+1)
      IF (ARG_PRESENT(history) EQ TRUE) THEN $
        history_spc = fltarr(n_candidate+1)
      IF (ARG_PRESENT(swept) EQ TRUE) THEN $
        swept_spc = fltarr(n_candidate+1)
   END
   ;
   ; Call the system function.
   ;
   err_status = 0L
   MATHSTAT_192, type, err_status, x_cvt, y_cvt, n_rows, n_candidate, $
                 all_steps_cvt, $
                 first_step_cvt, $
                 inter_step_cvt, $
                 last_step_cvt, $
                 forward_cvt, $
                 backward_cvt, $
                 stepwise_cvt, $
                 n_steps_cvt, $
                 force_cvt, $
                 cov_nobs_cvt, $
                 level_cvt, $
                 p_in_cvt, $
                 p_out_cvt, $
                 tolerance_cvt, $
                 frequencies_cvt, $
                 weights_cvt, $
                 cov_input_cvt, $
                 iend_spc, $
                 anova_table_spc, $
                 coef_t_spc, $
                 coef_vif_spc, $
                 cov_swept_spc, $
                 history_spc, $
                 swept_spc
   ;
   ; Now copy over all output keywords results.
   ;
      IF (ARG_PRESENT(anova_table) EQ TRUE) THEN $
        anova_table = anova_table_spc(0:12)
      IF (ARG_PRESENT(cov_swept) EQ TRUE) THEN $
        cov_swept=transpose(cov_swept_spc)  ;NOTE transpose()
      IF (ARG_PRESENT(coef_t_tests) EQ TRUE) THEN $
        coef_t_tests=transpose(coef_t_spc)  ;NOTE transpose()
      IF (ARG_PRESENT(coef_vif) EQ TRUE) THEN $
        coef_vif=coef_vif_spc 
      IF (ARG_PRESENT(history) EQ TRUE) THEN $
        history=history_spc
      IF (ARG_PRESENT(iend) EQ TRUE) THEN $
        iend = iend_spc
      IF (ARG_PRESENT(swept) EQ TRUE) THEN $
        swept=swept_spc
   ;
   ; Return.
   ;
   RETURN
END
