; $Id: //depot/idl/IDL_70/idl_src/libs/imsl/imsl_6.0/lib/imsl_robust_cov.pro#1 $
;
; Copyright (c) 1970-2006, VISUAL NUMERICS Inc. All Rights Reserved.
; This software is confidential information which is proprietary to and a
; trade secret of Visual Numerics Inc.  Use, duplication or disclosure is
; subject to the terms of an appropriate license agreement.
;
FUNCTION imsl_robust_cov, x, $                ;INPUT 2-D array: floating point
                     n_groups, $              ;INPUT Scalar LONG
                     init_est_mean=init_est_mean, $      ;INPUT Scalar ON/OFF flag
                     init_est_median=init_est_median, $  ;INPUT Scalar ON/OFF flag
                     mean_est=mean_est, $     ;INPUT 2-D Scalar floating point
                     cov_est=cov_est, $       ;INPUT 2-D Scalar floating point
                     idx_vars=idx_vars, $     ;INPUT 1-D array: LONG
                     idx_cols=idx_cols, $     ;INPUT 1-D array: LONG
                     stahel=stahel, $         ;INPUT Scalar ON/OFF flag
                     huber=huber, $           ;INPUT Scalar ON/OFF flag
                     percentage=percentage, $ ;INPUT Scalar floating point
                     itmax=itmax, $           ;INPUT Scalar LONG
                     tolerance=tolerance, $   ;INPUT Scalar floating point
                     minimax_weights=minimax_weights, $ ;OUTPUT 1-D array: floating point
                     group_counts=group_counts, $ ;OUTPUT 1-D array: LONG
                     sum_weights=sum_weights, $ ;OUTPUT 1-D array: floating point
                     means=means, $           ;OUTPUT 2-D array: floating point
                     u=u, $                   ;OUTPUT 2-D array: floating point
                     beta=beta, $             ;OUTPUT Scalar floating point                     
                     double=double, $         ;INPUT Scalar ON/OFF flag
                     nmissing=nmissing        ;OUTPUT Scalar LONG

@imsl_init.pro
   ON_ERROR, on_err_action

   ;
   ; Error checking:
   ; The following checks are performed.
   ;   X must be a 2D array (n_rows by n_var+1)
   ;   If supplied, IDX_VARS must be a 1D array of length nvars.
   ;                set nvars based on IDX_VARS, and XCOL_DIM based on X.
   ;   If supplied, IDX_COLS must be a 1D array of length 3.
   ;   Keywords MEAN_EST and COV_EST must be used together.
   ;   If supplied, MEAN_EST must be a 1D or 2D array, (n_groups by n_var)
   ;   If supplied, COV_EST must be a 1D or 2D array, (n_var by n_var)
   ;   Keywords MEAN_EST, INIT_EST_MEAN, and INIT_EST_MEDIAN are mutually exclusive.
   ;   Keywords STAHEL and HUBER are mutually exclusive.
   ;
   nargs = n_params()
   IF (nargs NE 2) THEN message, "Incorrect number of arguments."
   size_x = IMSL_SIZE(x)
   IF (size_x(0) NE 2) THEN BEGIN
      message, "X must be a 2-D array."
   END
   IF (size_x(2) LE  1) THEN BEGIN
      message, "The second dimension of X must be two or more."
   END
   n_groups_cvt  =  IMSL_LONG(n_groups(0))
   IF (n_groups_cvt LT 1) then MESSAGE, "N_GROUPS must be at least 1."
   n_rows = size_x(1)
   n_var   =   size_x(2)-1
   x_col_dim = size_x(2)
   ; If IDX_* keywords are supplied, then n_var is dependent on them.
   IF ((KEYWORD_SET(IDX_COLS) + KEYWORD_SET(IDX_VARS)) EQ 1) THEN $
     MESSAGE,  "Keywords IDX_VARS and IDX_COLS must be used together."
   IF (KEYWORD_SET(IDX_COLS)) THEN BEGIN
      size_idx_cols = IMSL_SIZE(idx_cols)
      IF (size_idx_cols(0) NE 1) THEN BEGIN
         message, "IDX_COLS must be a 1-D array."
      END
      IF (size_idx_cols(1) NE 3) THEN BEGIN
         message, "IDX_COLS is not the correct size."
      END
      size_idx_vars = IMSL_SIZE(idx_vars)
      IF (size_idx_vars(0) NE 1) THEN BEGIN
         message, "IDX_VARS must be a 1-D array."
      END
      n_var = IMSL_N_ELEMENTS(idx_vars)
   END
   IF ((KEYWORD_SET(mean_est) + KEYWORD_SET(cov_est)) EQ 1) THEN $
     MESSAGE,  "Keywords MEAN_SET and COV_EST must be used together."
   IF KEYWORD_SET(mean_est) THEN BEGIN
      size_mean_est = IMSL_SIZE(mean_est)
      IF (n_var GT 1) THEN BEGIN ; MEAN_EST must be 2D (n_groups x n_var)
         IF (size_mean_est(0) NE 2) THEN BEGIN
           message, "MEAN_EST must be a 2-D array."
         END
         IF ((size_mean_est(1) NE n_groups_cvt) OR (size_mean_est(2) NE n_var)) THEN $
           MESSAGE,   "MEAN_EST is not the correct size."
      END ELSE BEGIN ; MEAN_EST must be 1D (n_groups)
         IF (size_mean_est(0) NE 1) THEN BEGIN
           message, "MEAN_EST is not the correct size."
         END
         IF (size_mean_est(1) NE n_groups_cvt) THEN $
           MESSAGE,   "MEAN_EST is not the correct size."
      END
      size_cov_est = IMSL_SIZE(cov_est)
      IF (n_var GT 1) THEN BEGIN ; COV_EST must be 2D (n_var x n_var)
         IF (size_cov_est(0) NE 2) THEN BEGIN
           message, "COV_EST must be a 2-D array."
         END
         IF ((size_cov_est(1) NE n_var) OR (size_cov_est(2) NE n_var)) THEN $
           MESSAGE,   "COV_EST is not the correct size."
      END ELSE BEGIN ; COV_EST must be 1D (n_var)
         IF (size_cov_est(0) NE 1) THEN BEGIN
           message, "COV_EST is not the correct size."
         END
         IF (size_cov_est(1) NE n_var) THEN $
           MESSAGE,   "COV_EST is not the correct size."
      END
   END
   IF ((KEYWORD_SET(init_est_mean)+KEYWORD_SET(init_est_medean) + $
        KEYWORD_SET(mean_est)) GT 1) THEN $
      MESSAGE,   "Keywords INIT_EST_MEAN, INIT_EST_MEAN and MEAN_EST are mutually exclusive."
   IF ((KEYWORD_SET(stahel)+KEYWORD_SET(huber)) GT 1) THEN $
      MESSAGE,   "Keywords STAHEL and HUGBER are mutually exclusive."
   ;
   ; Decide on what precision to use.
   ; Use the highest precision of the input argument(s).
   ;
   type = TYP_FLOAT
   IF (size_x(N_ELEMENTS(size_x)-2) GT TYP_FLOAT) THEN type = TYP_DOUBLE
   IF (KEYWORD_SET(double) EQ true) THEN type = TYP_DOUBLE
   ;
   ; Setup the parameters for the call to the system function.
   ;
   ; Input LONG keyword(s)
   IF (KEYWORD_SET(itmax) EQ TRUE) THEN itmax_cvt  =  IMSL_LONG(ITMAX(0))
   IF (KEYWORD_SET(idx_vars) EQ TRUE) THEN idx_vars_cvt  =  IMSL_LONG(idx_vars)
   IF (KEYWORD_SET(idx_cols) EQ TRUE) THEN idx_cols_cvt  =  IMSL_LONG(idx_cols)
   IF (KEYWORD_SET(stahel) EQ TRUE) THEN stahel_cvt  =  IMSL_1 
   IF (KEYWORD_SET(huber) EQ TRUE) THEN huber_cvt  =  IMSL_1
   IF (KEYWORD_SET(init_est_mean) EQ TRUE) THEN i_e_mean_cvt  =  IMSL_1
   IF (KEYWORD_SET(init_est_median) EQ TRUE) THEN i_e_median_cvt  =  IMSL_1
   
   ; Output LONG keyword(s)
   IF (ARG_PRESENT(group_counts) EQ TRUE) THEN group_counts_spc = IMSL_LONARR(n_groups)
   IF (ARG_PRESENT(nmissing) EQ TRUE) THEN nmissing_spc = IMSL_LONG(0)
   ;
   ; Floating point arguments and keywords
   IF (type EQ TYP_DOUBLE) THEN BEGIN
      ; Input
      x_cvt = double(TRANSPOSE(x))
      IF (KEYWORD_SET(mean_est) EQ TRUE) THEN mean_est_cvt = DOUBLE(TRANSPOSE(mean_est))
      IF (KEYWORD_SET(cov_est) EQ TRUE) THEN cov_est_cvt = DOUBLE(TRANSPOSE(cov_est))
      IF (KEYWORD_SET(percentage) EQ TRUE) THEN percentage_cvt = DOUBLE(percentage(0)) ELSE percentage_cvt = DOUBLE(5.)
      IF (KEYWORD_SET(tolerance) EQ TRUE) THEN tolerance_cvt = DOUBLE(tolerance(0)) ELSE tolerance_cvt = DOUBLE(0.0001)
      ; Output
      IF (ARG_PRESENT(minimax_weights) EQ TRUE) THEN minmax_wts_spc  =  DBLARR(3)
      IF (ARG_PRESENT(beta) EQ TRUE) THEN beta_spc  =  DOUBLE(0.0)
      IF (ARG_PRESENT(sum_weights) EQ TRUE) THEN sum_weights_spc  =  DBLARR(n_groups_cvt)
      IF (ARG_PRESENT(means) EQ TRUE) THEN means_spc  =  DBLARR(n_var, n_groups)
      IF (ARG_PRESENT(u) EQ TRUE) THEN u_spc  =  DBLARR(n_var, n_var)
      result = dblarr(n_var, n_var)
   END ELSE BEGIN
      ; Input
      x_cvt = float(TRANSPOSE(x))
      IF (KEYWORD_SET(mean_est) EQ TRUE) THEN mean_est_cvt = FLOAT(TRANSPOSE(mean_est))
      IF (KEYWORD_SET(cov_est) EQ TRUE) THEN cov_est_cvt = FLOAT(TRANSPOSE(cov_est))
      IF (KEYWORD_SET(percentage) EQ TRUE) THEN percentage_cvt = FLOAT(percentage(0)) ELSE percentage_cvt = FLOAT(5.)
      IF (KEYWORD_SET(tolerance) EQ TRUE) THEN tolerance_cvt = FLOAT(tolerance(0)) ELSE tolerance_cvt = FLOAT(0.0001)
      ; Output
      IF (ARG_PRESENT(minimax_weights) EQ TRUE) THEN minmax_wts_spc  =  FLTARR(3)
      IF (ARG_PRESENT(beta) EQ TRUE) THEN beta_spc  =  FLOAT(0.0)
      IF (ARG_PRESENT(sum_weights) EQ TRUE) THEN sum_weights_spc  =  FLTARR(n_groups_cvt)
      IF (ARG_PRESENT(means) EQ TRUE) THEN means_spc  =  FLTARR(n_var, n_groups)
      IF (ARG_PRESENT(u) EQ TRUE) THEN u_spc  =  FLTARR(n_var, n_var)
      result = fltarr(n_var, n_var)
   END
   ;
   ; Call the system function.
   ;
   err_status = 0L

   MATHSTAT_280,type,err_status, $
                           x_cvt, $
                           x_col_dim, $
                           n_var, $
                           n_rows, $
                           n_groups_cvt, $
                           itmax_cvt, $
                           idx_vars_cvt, $
                           idx_cols_cvt, $
                           stahel_cvt, $
                           huber_cvt, $
                           i_e_mean_cvt, $
                           i_e_median_cvt, $
                           percentage_cvt, $
                           tolerance_cvt, $
                           mean_est_cvt, $
                           cov_est_cvt, $
                           minmax_wts_spc, $
                           beta_spc, $
                           sum_weights_spc, $
                           means_spc, $
                           u_spc, $
                           group_counts_spc, $
                           nmissing_spc, $
                           result
   ;
   ; Now copy over all output keywords results.
   ;
      IF (ARG_PRESENT(minimax_weights) EQ TRUE) THEN minimax_weights = minmax_wts_spc
      IF (ARG_PRESENT(beta) EQ TRUE) THEN beta = beta_spc
      IF (ARG_PRESENT(sum_weights) EQ TRUE) THEN sum_weights = sum_weights_spc
      IF (ARG_PRESENT(means) EQ TRUE) THEN means = TRANSPOSE(means_spc)
      IF (ARG_PRESENT(u) EQ TRUE) THEN u = TRANSPOSE(u_spc)
      IF (ARG_PRESENT(group_counts) EQ TRUE) THEN group_counts = group_counts_spc
      IF (ARG_PRESENT(nmissing) EQ TRUE) THEN nmissing = nmissing_spc
   ;
   ; Return.
   ;
   RETURN, TRANSPOSE(result)
END

                   
                   
                   

  
      

  
