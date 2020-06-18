; $Id: //depot/idl/IDL_70/idl_src/libs/imsl/imsl_6.0/lib/imsl_factor_analysis.pro#1 $
;
; Copyright (c) 1970-2006, VISUAL NUMERICS Inc. All Rights Reserved.
; This software is confidential information which is proprietary to and a
; trade secret of Visual Numerics Inc.  Use, duplication or disclosure is
; subject to the terms of an appropriate license agreement.
;
FUNCTION imsl_factor_analysis, covariances, $              ;INPUT 1-D array: floating point
                          n_factors, $                     ;INPUT Scalar LONG
                          Max_Likelihood=m_ml, $           ;INPUT Scalar ON/OFF flag
                          Princ_Comp=m_cp, $               ;INPUT Scalar ON/OFF flag
                          Princ_Factor=m_pf, $             ;INPUT Scalar ON/OFF flag
                          Unwgt_Lsq=m_uls, $               ;INPUT Scalar ON/OFF flag
                          Gen_Lsq=m_ls, $                  ;INPUT Scalar ON/OFF flag
                          Image=m_image, $                 ;INPUT Scalar ON/OFF flag
                          Alpha=m_alpha, $                 ;INPUT Scalar LONG
                          Unique_Var_In=unique_var_in, $   ;INPUT 1-D array: floating point
                          Itmax=itmax, $                   ;INPUT Scalar LONG
                          Max_Steps=max_steps, $           ;INPUT Scalar LONG
                          Eps=eps, $                       ;INPUT Scalar floating point
                          Switch_Eps=switch_eps, $         ;INPUT Scalar floating point
                          Double=double, $                 ;INPUT Scalar ON/OFF flag
                          Iters=iters, $                   ;OUTPUT Scalar LONG
                          Unique_Var_out=unique_var_out, $ ;OUTPUT 1-D array: floating point
                          Eigenvalues=eigenvalues, $       ;OUTPUT 1-D array: floating point
                          Chi_Sq_Test=chi_sq_test, $       ;OUTPUT 1-D array: floating point
                          Last_Step=last_step, $           ;OUTPUT 1-D array: floating point
                          Tucker_Coef=tucker_coef, $       ;OUTPUT Scalar floating point
                          F_Min=f_min                      ;OUTPUT Scalar floating point
  
@imsl_init.pro
   ON_ERROR, on_err_action
   ;
   ; Error checking.
   ; Make sure the input arguments are 1-D or 2-D arrays.
   ;
   nargs = n_params()
   IF (nargs NE 2) THEN $
         message, "Incorrect number of arguments."
   size_c = IMSL_SIZE(covariances)
   IF (size_c(0) NE 2) THEN $
      message, "The input array, COVARIANCES, must be 2-D"
   IF (size_c(1) NE size_c(2)) THEN $
      message, "The input array, COVARIANCES, must be square"
   nvar = size_c(1)
   IF (KEYWORD_SET(unique_var_in)) THEN BEGIN
      size_unique_var_in = IMSL_SIZE(unique_var_in)
      IF (size_unique_var_in(0) NE 1) THEN $
         message, "The input array, UNIQUE_VAR_IN, must be 1-D"
      IF (size_unique_var_in(1) NE nvar) THEN $ 
         message, "The length of the UNIQUE_VAR_IN array  be equal to the number of variables."
   END
   IF (n_factors LT 1) THEN $
      message, "The number of factors, N_FACTORS, must be positive."
   ;
   ; Decide on what precision to use.
   ; Use the highest precision of the input argument(s).
   ;
   type = TYP_FLOAT
   IF (size_c(N_ELEMENTS(size_c)-2) GT TYP_FLOAT) THEN type = TYP_DOUBLE
   IF (KEYWORD_SET(double) EQ true) THEN type = TYP_DOUBLE
   ;
   ; Setup the parameters for the call to the system function.
   ;
   ; Input LONG argument(s)
   n_factors_cvt = (IMSL_LONG(n_factors))(0)
   ;
   ; Input LONG keyword(s)
   IF (KEYWORD_SET(m_ml) EQ TRUE) THEN $
     m_ml_cvt = (IMSL_LONG(m_ml))(0)
   IF (KEYWORD_SET(m_pc) EQ TRUE) THEN $
     m_pc_cvt = IMSL_1
   IF (KEYWORD_SET(m_pf) EQ TRUE) THEN $
     m_pf_cvt = IMSL_1
   IF (KEYWORD_SET(m_uls) EQ TRUE) THEN $
     m_uls_cvt = IMSL_1
   IF (KEYWORD_SET(m_ls) EQ TRUE) THEN $
     m_ls_cvt = (IMSL_LONG(m_ls))(0)
   IF (KEYWORD_SET(m_image) EQ TRUE) THEN $
     m_image_cvt = IMSL_1
   IF (KEYWORD_SET(m_alpha) EQ TRUE) THEN $
     m_alpha_cvt = (IMSL_LONG(m_alpha))(0)
   IF (KEYWORD_SET(itmax) EQ TRUE) THEN $
     itmax_cvt = (IMSL_LONG(itmax))(0)
   IF (KEYWORD_SET(max_steps) EQ TRUE) THEN $
     max_steps_cvt = (IMSL_LONG(max_steps))(0)
   ;
   ; Output LONG keyword(s)
   IF (ARG_PRESENT(iters) EQ TRUE) THEN $
     iters_spc = IMSL_0
   ;
   ; Floating point arguments and keywords
   IF (type EQ TYP_DOUBLE) THEN BEGIN
      ; Result vector.
      l_result = dblarr(n_factors, nvar)
      ; 
      ; Input 
      covariances_cvt = double(covariances)
      IF (KEYWORD_SET(unique_var_in) EQ TRUE) THEN $
        u_var_in_cvt = double(unique_var_in)
      IF (KEYWORD_SET(eps) EQ TRUE) THEN $
        eps_cvt = (double(eps))(0)
      IF (KEYWORD_SET(switch_eps) EQ TRUE) THEN $
        switch_eps_cvt = (double(switch_eps))(0)
      ;
      ; Output keywords.
      ;
      IF (ARG_PRESENT(unique_var_out) EQ TRUE) THEN $
        u_var_out_spc = dblarr(nvar)
      IF (ARG_PRESENT(eigenvalues) EQ TRUE) THEN $
        eigenvalues_spc = dblarr(nvar)
      IF (ARG_PRESENT(chi_sq_test) EQ TRUE) THEN $
        chi_sq_test_spc = dblarr(3)
      IF (ARG_PRESENT(tucker_coef) EQ TRUE) THEN $
        tucker_coef_spc = double(0.0)
      IF (ARG_PRESENT(f_min) EQ TRUE) THEN $
        f_min_spc = double(0.0)
      IF (ARG_PRESENT(last_step) EQ TRUE) THEN $
        last_step_spc = dblarr(nvar)
   END ELSE BEGIN
      ; Result vector.
      l_result = fltarr(n_factors, nvar)
      ; 
      ; Input 
      covariances_cvt = float(covariances)
      IF (KEYWORD_SET(unique_var_in) EQ TRUE) THEN $
        u_var_in_cvt = float(unique_var_in)
      IF (KEYWORD_SET(eps) EQ TRUE) THEN $
        eps_cvt = (float(eps))(0)
      IF (KEYWORD_SET(switch_eps) EQ TRUE) THEN $
        switch_eps_cvt = (float(switch_eps))(0)
      ;
      ; Output keywords.
      ;
      IF (ARG_PRESENT(unique_var_out) EQ TRUE) THEN $
        u_var_out_spc = fltarr(nvar)
      IF (ARG_PRESENT(eigenvalues) EQ TRUE) THEN $
        eigenvalues_spc = fltarr(nvar)
      IF (ARG_PRESENT(chi_sq_test) EQ TRUE) THEN $
        chi_sq_test_spc = fltarr(3)
      IF (ARG_PRESENT(tucker_coef) EQ TRUE) THEN $
        tucker_coef_spc = float(0.0)
      IF (ARG_PRESENT(f_min) EQ TRUE) THEN $
        f_min_spc = float(0.0)
      IF (ARG_PRESENT(last_step) EQ TRUE) THEN $
        last_step_spc = fltarr(nvar)
   END
   ;
   ; Call the system function.
   ;
   err_status = 0L
   MATHSTAT_137, type, err_status, covariances_cvt, nvar, n_factors_cvt, $
     m_ml_cvt, m_pc_cvt, m_pf_cvt, m_uls_cvt, m_ls_cvt, m_image_cvt, m_alpha_cvt, $
     u_var_in_cvt, itmax_cvt, max_steps_cvt, eps_cvt, switch_eps_cvt,  $
     u_var_out_spc, eigenvalues_spc, chi_sq_test_spc, $
     tucker_coef_spc, iters_spc, f_min_spc, last_step_spc, l_result 
   ;
   ; Now copy over all output keywords results.
   ;
   IF (ARG_PRESENT(chi_sq_test) EQ TRUE) THEN $
     chi_sq_test=chi_sq_test_spc
   IF (ARG_PRESENT(eigenvalues) EQ TRUE) THEN $
     eigenvalues=eigenvalues_spc
   IF (ARG_PRESENT(f_min) EQ TRUE) THEN $
     f_min=f_min_spc
   IF (ARG_PRESENT(iters) EQ TRUE) THEN $
     iters=iters_spc
   IF (ARG_PRESENT(last_step) EQ TRUE) THEN $
     last_step=last_step_spc
   IF (ARG_PRESENT(unique_var_out) EQ TRUE) THEN $
     unique_var_out=u_var_out_spc
   IF (ARG_PRESENT(tucker_coef) EQ TRUE) THEN $
     tucker_coef=tucker_coef_spc
   ;
   ; Return.
   ;
   RETURN, transpose(l_result)
END
   

