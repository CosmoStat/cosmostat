; $Id: //depot/idl/IDL_70/idl_src/libs/imsl/imsl_6.0/lib/imsl_discr_analysis.pro#1 $
;
; Copyright (c) 1970-2006, VISUAL NUMERICS Inc. All Rights Reserved.
; This software is confidential information which is proprietary to and a
; trade secret of Visual Numerics Inc.  Use, duplication or disclosure is
; subject to the terms of an appropriate license agreement.
;
PRO imsl_discr_analysis, x, $                   ;INPUT 2-D array: floating point
                     n_groups, $                ;INPUT Scalar LONG
                     idx_vars=idx_vars, $       ;INPUT 1-D array: LONG
                     idx_cols=idx_cols, $       ;INPUT 1-D array: LONG
                     method=method, $           ;INPUT Scalar LONG
                     prior_equal=prior_equal, $ ;INPUT Scalar ON/OFF flag
                     prior_prop=prior_prop, $   ;INPUT Scalar ON/OFF flag
                     prior_input=prior_input, $ ;INPUT 1-D array: floating point
                     prior_output=prior_output, $ ;OUTPUT 1-D array: floating point
                     group_counts=group_counts, $ ;OUTPUT 1-D array: LONG
                     means=means, $               ;OUTPUT 2-D array: floating point
                     covariances=covariances, $   ;OUTPUT 3-D array: floating point
                     coefficients=coefficients, $ ;OUTPUT 2-D array: floating point
                     class_member=class_member, $ ;OUTPUT 1-D array: LONG
                     class_table=class_table, $   ;OUTPUT 2-D array: floating point
                     prob=prob, $                 ;OUTPUT 2-D array: floating point
                     mahalanobis=mahalanobis, $   ;OUTPUT 2-D array: floating point
                     stats=stats, $               ;OUTPUT 1-D array: floating point
                     double=double, $             ;INPUT Scalar ON/OFF flag
                     nmissing=nmissing            ;OUTPUT Scalar LONG

@imsl_init.pro
   ON_ERROR, on_err_action

   ;
   ; Error checking:
   ; The following checks are performed.
   ;   X must be a 2D array (n_rows by n_var+1)
   ;   If supplied, IDX_VARS must be a 1D array of length n_var
   ;                set n_var based on IDX_VARS, and XCOL_DIM based on X.
   ;   If supplied, IDX_COLS must be a 1D array of length 3.
   ;   Keywords PRIOR_* are mutually exclusive.
   ;   If supplied, PRIOR_INPUT must be a 1D array, (n_groups)
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
   IF ((KEYWORD_SET(prior_equal)+KEYWORD_SET(prior_prop) + $
        KEYWORD_SET(prior_input)) GT 1) THEN $
      MESSAGE,   "Keywords PRIOR_EQUAL, PRIOR_PROP and PRIOR_INPUT are mutually exclusive."
   IF (KEYWORD_SET(prior_input)) THEN BEGIN
      size_prior_input = IMSL_SIZE(prior_input)
      IF (size_prior_input(0) NE 1) THEN BEGIN
         message, "PRIOR_INPUT must be a 1-D array."
      END
      IF (size_prior_input(1) NE n_groups_cvt) THEN $
        MESSAGE,   "Keyword PRIOR_INPUT is not the correct size."
   END
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
   IF (KEYWORD_SET(idx_vars) EQ TRUE) THEN idx_vars_cvt  =  IMSL_LONG(idx_vars)
   IF (KEYWORD_SET(idx_cols) EQ TRUE) THEN idx_cols_cvt  =  IMSL_LONG(idx_cols)
   IF (KEYWORD_SET(method) EQ TRUE) THEN method_cvt  =  IMSL_LONG(method(0)) ELSE method_cvt = IMSL_1 
   IF (KEYWORD_SET(prior_equal) EQ TRUE) THEN prior_equal_cvt  =  IMSL_1
   IF (KEYWORD_SET(prior_prop) EQ TRUE) THEN prior_prop_cvt  =  IMSL_1
   
   ; Output LONG keyword(s)
   IF (ARG_PRESENT(group_counts) EQ TRUE) THEN group_counts_spc = IMSL_LONARR(n_groups)
   IF (ARG_PRESENT(class_member) EQ TRUE) THEN class_member_spc = IMSL_LONARR(n_rows)
   IF (ARG_PRESENT(nmissing) EQ TRUE) THEN nmissing_spc = IMSL_LONG(0)
   ;
   IF ((method_cvt EQ 3)  OR (method_cvt EQ 6)) THEN g = 1 ELSE g = n_groups_cvt+1
   ; Floating point arguments and keywords
   IF (type EQ TYP_DOUBLE) THEN BEGIN
      ; Input
      x_cvt = double(TRANSPOSE(x))
      IF (KEYWORD_SET(prior_input) EQ TRUE) THEN prior_input_cvt = DOUBLE(prior_input)
      ; Output
      IF (ARG_PRESENT(prior_output) EQ TRUE) THEN prior_output_spc  =  DBLARR(n_groups_cvt)
      IF (ARG_PRESENT(means) EQ TRUE) THEN means_spc  =  DBLARR(n_var, n_groups_cvt)
      IF (ARG_PRESENT(covariances) EQ TRUE) THEN covariances_spc  =  DBLARR(n_var, n_var, g)
      IF (ARG_PRESENT(coefficients) EQ TRUE) THEN coefficients_spc  =  DBLARR(n_var+1, n_groups_cvt)
      IF (ARG_PRESENT(class_table) EQ TRUE) THEN class_table_spc  =  DBLARR(n_groups_cvt, n_groups_cvt)
      IF (ARG_PRESENT(prob) EQ TRUE) THEN prob_spc  =  DBLARR(n_groups_cvt, n_rows)
      IF (ARG_PRESENT(mahalanobis) EQ TRUE) THEN mahalanobis_spc  =  DBLARR(n_groups_cvt, n_groups_cvt)
      IF (ARG_PRESENT(stats) EQ TRUE) THEN stats_spc  =  DBLARR(IMSL_4+IMSL_2*(n_groups_cvt+1))
   END ELSE BEGIN
      ; Input
      x_cvt = float(TRANSPOSE(x))
      IF (KEYWORD_SET(prior_input) EQ TRUE) THEN prior_input_cvt = FLOAT(prior_input)
      ; Output
      IF (ARG_PRESENT(prior_output) EQ TRUE) THEN prior_output_spc  =  FLTARR(n_groups_cvt)
      IF (ARG_PRESENT(means) EQ TRUE) THEN means_spc  =  FLTARR(n_var, n_groups_cvt)
      IF (ARG_PRESENT(covariances) EQ TRUE) THEN covariances_spc  =  FLTARR(n_var, n_var, g)
      IF (ARG_PRESENT(coefficients) EQ TRUE) THEN coefficients_spc  =  FLTARR(n_var+1, n_groups_cvt)
      IF (ARG_PRESENT(class_table) EQ TRUE) THEN class_table_spc  =  FLTARR(n_groups_cvt, n_groups_cvt)
      IF (ARG_PRESENT(prob) EQ TRUE) THEN prob_spc  =  FLTARR(n_groups_cvt, n_rows)
      IF (ARG_PRESENT(mahalanobis) EQ TRUE) THEN mahalanobis_spc  =  FLTARR(n_groups_cvt, n_groups_cvt)
      IF (ARG_PRESENT(stats) EQ TRUE) THEN stats_spc  =  FLTARR(IMSL_4+IMSL_2*(n_groups_cvt+1))
   END
   ;
   ; Call the system function.
   ;
   err_status = 0L
   MATHSTAT_281,type,err_status, $
                           x_cvt, $
                           x_col_dim, $
                           n_var, $
                           n_rows, $
                           n_groups_cvt, $
                           idx_vars_cvt, $
                           idx_cols_cvt, $
                           method_cvt, $
                           prior_equal_cvt, $
                           prior_prop_cvt, $
                           prior_input_cvt, $
                           prior_output_spc, $
                           group_counts_spc, $
                           means_spc, $
                           covariances_spc, $
                           coefficients_spc, $
                           class_member_spc, $
                           class_table_spc, $
                           prob_spc, $
                           mahalanobis_spc, $
                           stats_spc, $
                           nmissing_spc
   ;
   ; Now copy over all output keywords results.
   ;
   IF (ARG_PRESENT(prior_output) EQ TRUE) THEN prior_output = prior_output_spc
   IF (ARG_PRESENT(group_counts) EQ TRUE) THEN group_counts = group_counts_spc
   IF (ARG_PRESENT(class_member) EQ TRUE) THEN class_member = class_member_spc
   IF (ARG_PRESENT(nmissing) EQ TRUE) THEN nmissing = nmissing_spc
   IF (ARG_PRESENT(means) EQ TRUE) THEN means = TRANSPOSE(means_spc)
   IF (ARG_PRESENT(coefficients) EQ TRUE) THEN coefficients = TRANSPOSE(coefficients_spc)
   IF (ARG_PRESENT(class_table) EQ TRUE) THEN class_table = TRANSPOSE(class_table_spc)
   IF (ARG_PRESENT(prob) EQ TRUE) THEN prob = TRANSPOSE(prob_spc)
   IF (ARG_PRESENT(mahalanobis) EQ TRUE) THEN mahalanobis = TRANSPOSE(mahalanobis_spc)
   IF (ARG_PRESENT(stats) EQ TRUE) THEN stats = stats_spc
   IF (ARG_PRESENT(covariances) EQ TRUE) THEN BEGIN
      IF (type EQ TYP_DOUBLE) THEN covariances = DBLARR(g, n_var,n_var) $
        ELSE  covariances = FLTARR(g, n_var,n_var) 
      FOR i = 0, n_var-1 DO FOR j   =   0,   n_var-1 DO $
         covariances(*,  j,  i) = covariances_spc(i, j, *)
   END
   ;
   ; Return.
   ;
   RETURN
END

                   
                   
                   

  
      

  
