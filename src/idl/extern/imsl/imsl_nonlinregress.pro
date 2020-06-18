; $Id: //depot/idl/IDL_70/idl_src/libs/imsl/imsl_6.0/lib/imsl_nonlinregress.pro#1 $
;
; Copyright (c) 1970-2006, VISUAL NUMERICS Inc. All Rights Reserved.
; This software is confidential information which is proprietary to and a
; trade secret of Visual Numerics Inc.  Use, duplication or disclosure is
; subject to the terms of an appropriate license agreement.
;
FUNCTION imsl_nonlinregress, fcn, $               ;INPUT 1-D array: floating point
               n_parameters, $                    ;INPUT Scalar LONG
               x, $                               ;INPUT 1-D array: floating point
               y, $                               ;INPUT 1-D array: floating point
               itmax=itmax, $                     ;INPUT Scalar LONG
               n_digit=n_digit, $                 ;INPUT Scalar LONG
               max_jac_evals=max_jac_evals, $     ;INPUT Scalar LONG
               max_sse_evals=max_sse_evals, $     ;INPUT Scalar LONG
               double=double, $                   ;INPUT Scalar ON/OFF flag
               jacobian=jacobian, $               ;INPUT Scalar STRING
               grad_eps=grad_eps, $               ;INPUT Scalar floating point
               step_eps=step_eps, $               ;INPUT Scalar floating point
               rel_eps_sse=rel_eps_sse, $         ;INPUT Scalar floating point
               abs_eps_sse=abs_eps_sse, $         ;INPUT Scalar floating point
               max_step=max_step, $               ;INPUT Scalar floating point
               trust_region=trust_region, $       ;INPUT Scalar floating point
               tolerance=tolerance, $             ;INPUT Scalar floating point
               theta_guess=theta_guess, $         ;INPUT 1-D array: floating point
               theta_scale=theta_scale, $         ;INPUT 1-D array: floating point
               df=df, $                           ;OUTPUT Scalar LONG
               r_rank=r_rank, $                   ;OUTPUT Scalar LONG
               sse=sse, $                         ;OUTPUT Scalar floating point
               predicted=predicted, $             ;OUTPUT 1-D array: floating point
               residual=residual, $               ;OUTPUT 1-D array: floating point
               r_matrix=r_matrix                  ;OUTPUT 2-D array: floating point
               
@imsl_init.pro
   ON_ERROR, on_err_action
   ;
   ; Error checking.  
   ; The following checks are performed.
   ;   o  Checks on FCN
   ;      Must be a scalar string.
   ;   o  Checks on X
   ;      1. Must be an array.
   ;      2. Can be either a 2-D array of dimension (nobs)x(n_indep)
   ;         of a 1-D array of dimension (nobs).
   ;   o  Checks on p_argv[3]:
   ;      1. must be a 1-D array of length p_argv[2]->value.arr->dim[0].
   ;   o  If THETA_GUESS:
   ;      1. Must be an array of length n_params, if (n_params==1), THETA_GUESS can 
   ;         be a scalar.
   ;      2. Convert it to rtyp.
   ;   o  If THETA_SCALE:
   ;      1. Must be an array of length n_params, if (n_params==1), THETA_SCALE can 
   ;         be a scalar.
   ;
   nargs = n_params()
   IF (nargs NE 4) THEN $
         message, "Incorrect number of arguments."
   ;
   size_fcn = IMSL_SIZE(fcn)
   IF ((size_fcn(0) NE 0) OR (size_fcn(n_elements(size_fcn)-2) NE 7)) THEN BEGIN 
      message, "FCN must be a scalar string."
   END
   size_x = IMSL_SIZE(x)
   IF ((size_x(0) LT 1) or (size_x(0) GT 2)) THEN BEGIN 
      message, "The input array, X, must be 1-D or 2-D."
   END
   nobs = size_x(1) ;Define nobs
   IF (size_x(0) EQ 1) THEN n_indep = IMSL_1 ELSE n_indep = IMSL_LONG(size_x(2)) ; Define n_indep
   ;
   size_y = IMSL_SIZE(y)
   IF (size_y(0) NE 1)  THEN BEGIN 
      message, "The input array, Y, must be 1-D."
   END
   IF (size_y(1) NE nobs)  THEN BEGIN 
      message, "The input array, Y, is not the correct size."
   END
   ;
   n_params_cvt = (IMSL_LONG(n_parameters))(0)
   IF (KEYWORD_SET(theta_guess) EQ true) THEN BEGIN
      size_tmp = IMSL_SIZE(theta_guess)
      IF ((n_params_cvt GT 1) AND (size_tmp(0) NE 1)) THEN $
        message, "THETA_GUESS must be a 1-D array of length N_PARAMETERS."
      IF ((n_params_cvt GT 1) AND (size_tmp(1) NE n_params_cvt)) THEN $
        message, "THETA_GUESS must be a 1-D array of length N_PARAMETERS."      
      IF ((n_params_cvt EQ 1) AND (N_ELEMENTS(theta_guess) NE 1)) THEN $
        message, "THETA_GUESS must have exactly N_PARAMETERS elements."
   END
   n_params_cvt = (IMSL_LONG(n_parameters))(0)
   IF (KEYWORD_SET(theta_scale) EQ true) THEN BEGIN
      size_tmp = IMSL_SIZE(theta_scale)
      IF ((n_params_cvt GT 1) AND (size_tmp(0) NE 1)) THEN $
        message, "THETA_SCALE must be a 1-D array of length N_PARAMETERS."
      IF ((n_params_cvt GT 1) AND (size_tmp(1) NE n_params_cvt)) THEN $
        message, "THETA_SCALE must be a 1-D array of length N_PARAMETERS."      
      IF ((n_params_cvt EQ 1) AND (N_ELEMENTS(theta_scale) NE 1)) THEN $
        message, "THETA_SCALE must have exactly N_PARAMETERS elements."
   END
   ;
   ;ERROR CHECKING COMPLETE.
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
   IF (KEYWORD_SET(itmax) EQ TRUE) THEN $
     itmax_cvt = (IMSL_LONG(itmax))(0)
   IF (KEYWORD_SET(n_digit) EQ TRUE) THEN $
     n_digit_cvt = (IMSL_LONG(n_digit))(0)
   IF (KEYWORD_SET(max_jac_evals) EQ TRUE) THEN $
     max_jac_ev_cvt = (IMSL_LONG(max_jac_evals))(0)
   IF (KEYWORD_SET(max_sse_evals) EQ TRUE) THEN $
     max_sse_ev_cvt = (IMSL_LONG(max_sse_evals))(0)
   ;
   ; Input STRING keyword(s)
   IF (KEYWORD_SET(jacobian) EQ TRUE) THEN $
     jacobian_cvt = (STRING(jacobian))(0)
   ;
   ; Output LONG keyword(s)
   IF (ARG_PRESENT(df) EQ TRUE) THEN $
     df_spc = IMSL_0
   IF (ARG_PRESENT(r_rank) EQ TRUE) THEN $
     r_rank_spc = IMSL_0
   ;
   ; Floating point arguments and keywords
   IF (type EQ TYP_DOUBLE) THEN BEGIN
      ; Result vector.
      result = dblarr(n_params_cvt)
      ; 
      ; Input
      x_cvt = imsl_cvt_arr(x,type)
      y_cvt = double(y)
      IF (KEYWORD_SET(grad_eps) EQ TRUE) THEN $
        grad_eps_cvt = (DOUBLE(grad_eps))(0)
      IF (KEYWORD_SET(step_eps) EQ TRUE) THEN $
        step_eps_cvt = (DOUBLE(step_eps))(0)
      IF (KEYWORD_SET(rel_eps_sse) EQ TRUE) THEN $
        rel_eps_sse_cvt = (DOUBLE(rel_eps_sse))(0)
      IF (KEYWORD_SET(abs_eps_sse) EQ TRUE) THEN $
        abs_eps_sse_cvt = (DOUBLE(abs_eps_sse))(0)
      IF (KEYWORD_SET(max_step) EQ TRUE) THEN $
        max_step_cvt = (DOUBLE(max_step))(0)
      IF (KEYWORD_SET(trust_region) EQ TRUE) THEN $
        trst_region_cvt = (DOUBLE(trust_region))(0)
      IF (KEYWORD_SET(tolerance) EQ TRUE) THEN $
        tolerance_cvt = (DOUBLE(tolerance))(0)
      IF (KEYWORD_SET(theta_guess) EQ TRUE) THEN $
        theta_guess_cvt = (DOUBLE(theta_guess))
      IF (KEYWORD_SET(theta_scale) EQ TRUE) THEN $
        theta_scale_cvt = (DOUBLE(theta_scale))
      ; 
      ; Output 
      IF (ARG_PRESENT(sse) EQ TRUE) THEN $
         sse_spc = (DOUBLE(sse))(0)
      IF (ARG_PRESENT(predicted) EQ TRUE) THEN $
         predicted_spc = dblarr(nobs)
      IF (ARG_PRESENT(residual) EQ TRUE) THEN $
         residual_spc = dblarr(nobs)
      IF (ARG_PRESENT(r_matrix) EQ TRUE) THEN $
         r_matrix_spc = dblarr(n_params_cvt, n_params_cvt)
   END ELSE BEGIN
      ; Result vector.
      result = fltarr(n_params_cvt)
      ; 
      ; Input
      x_cvt = imsl_cvt_arr(x,type)
      y_cvt = float(y)
      IF (KEYWORD_SET(grad_eps) EQ TRUE) THEN $
        grad_eps_cvt = (FLOAT(grad_eps))(0)
      IF (KEYWORD_SET(step_eps) EQ TRUE) THEN $
        step_eps_cvt = (FLOAT(step_eps))(0)
      IF (KEYWORD_SET(rel_eps_sse) EQ TRUE) THEN $
        rel_eps_sse_cvt = (FLOAT(rel_eps_sse))(0)
      IF (KEYWORD_SET(abs_eps_sse) EQ TRUE) THEN $
        abs_eps_sse_cvt = (FLOAT(abs_eps_sse))(0)
      IF (KEYWORD_SET(max_step) EQ TRUE) THEN $
        max_step_cvt = (FLOAT(max_step))(0)
      IF (KEYWORD_SET(trust_region) EQ TRUE) THEN $
        trst_region_cvt = (FLOAT(trust_region))(0)
      IF (KEYWORD_SET(tolerance) EQ TRUE) THEN $
        tolerance_cvt = (FLOAT(tolerance))(0)
      IF (KEYWORD_SET(theta_guess) EQ TRUE) THEN $
        theta_guess_cvt = (FLOAT(theta_guess))
      IF (KEYWORD_SET(theta_scale) EQ TRUE) THEN $
        theta_scale_cvt = (FLOAT(theta_scale))
      ; 
      ; Output 
      IF (ARG_PRESENT(sse) EQ TRUE) THEN $
         sse_spc = (FLOAT(sse))(0)
      IF (ARG_PRESENT(predicted) EQ TRUE) THEN $
         predicted_spc = fltarr(nobs)
      IF (ARG_PRESENT(residual) EQ TRUE) THEN $
         residual_spc = fltarr(nobs)
      IF (ARG_PRESENT(r_matrix) EQ TRUE) THEN $
         r_matrix_spc = fltarr(n_params_cvt, n_params_cvt)
   END
   ;
   ; Call the system function.
   ;
   err_status = 0L
   MATHSTAT_165,   type, err_status, $
               fcn, $
               n_params_cvt, $
               x_cvt, $
               y_cvt, $
               nobs, $
               n_indep, $
               itmax_cvt, $
               n_digit_cvt, $
               max_jac_ev_cvt, $
               max_sse_ev_cvt, $
               jacobian_cvt, $
               grad_eps_cvt, $
               step_eps_cvt, $
               rel_eps_sse_cvt, $
               abs_eps_sse_cvt, $
               max_step_cvt, $
               trst_region_cvt, $
               tolerance_cvt, $
               theta_guess_cvt, $
               theta_scale_cvt, $
               df_spc, $
               r_rank_spc, $
               sse_spc, $
               predicted_spc, $
               residual_spc, $
               r_matrix_spc, $
               result
   ;
   ; Now copy over all output keywords results.
   ;
   IF (ARG_PRESENT(df) EQ TRUE) THEN $
     df = df_spc
   IF (ARG_PRESENT(r_rank) EQ TRUE) THEN $
     r_rank = r_rank_spc
   IF (ARG_PRESENT(sse) EQ TRUE) THEN $
     sse = sse_spc
   IF (ARG_PRESENT(predicted) EQ TRUE) THEN $
     predicted = predicted_spc
   IF (ARG_PRESENT(residual) EQ TRUE) THEN $
     residual = residual_spc
   IF (ARG_PRESENT(r_matrix) EQ TRUE) THEN $
     r_matrix = transpose(r_matrix_spc) ;NOTE transpose()
   ;
   ; Return.
   ;
   RETURN, result
END
   

