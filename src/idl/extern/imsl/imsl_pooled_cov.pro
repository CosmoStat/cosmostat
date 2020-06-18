; $Id: //depot/idl/IDL_70/idl_src/libs/imsl/imsl_6.0/lib/imsl_pooled_cov.pro#1 $
;
; Copyright (c) 1970-2006, VISUAL NUMERICS Inc. All Rights Reserved.
; This software is confidential information which is proprietary to and a
; trade secret of Visual Numerics Inc.  Use, duplication or disclosure is
; subject to the terms of an appropriate license agreement.
;
FUNCTION imsl_pooled_cov, X, $              ;INPUT 2-D array: floating point 
           ngroups,         $               ;INPUT Scalar LONG
           idx_vars=idx_vars, $             ;INPUT 1-D array: LONG
           idx_cols=idx_cols, $             ;INPUT 1-D array: LONG
           sum_weights=sum_weights, $       ;OUTPUT 1-D array: floating point 
           nmissing=nmissing, $             ;OUTPUT Scalar LONG
           gcounts=gcounts, $               ;OUTPUT 1-D array: floating point 
           means=means, $                   ;OUTPUT 2-D array: floating point 
           u=u, $                           ;OUTPUT 2-D array: floating point 
           double=double                    ;INPUT Scalar ON/OFF flag
           
@imsl_init.pro
   ON_ERROR, on_err_action
   ;
   ; Error checking:
   ;  - X must be a 2-D array. (m x n)
   ;  - NGROUPS must be positive.     
   nargs = n_params()
   IF (nargs NE 2)  THEN message, 'Incorrect number of arguments.'
     
   size_x = IMSL_LONG(size(x))
   IF ((size_x(0) NE 2)) THEN $
     message, 'X must be a 2-D array.'
   m = IMSL_LONG(size_x(1))
   n = IMSL_LONG(size_x(2)) - 1
   x_col_dim = size_x(2)
   IF (n LT 1) THEN MESSAGE, "The second dimension of X must be at least two."
   IF (ngroups LT 1) THEN MESSAGE, "NGROUPS must be positive."
   ; If IDX_* keywords are supplied, then n is dependent on them.
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
      n = IMSL_N_ELEMENTS(idx_vars)
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
   ng_cvt = IMSL_LONG(ngroups(0))
   IF (KEYWORD_SET(idx_vars) EQ TRUE) THEN idx_vars_cvt  =  IMSL_LONG(idx_vars)
   IF (KEYWORD_SET(idx_cols) EQ TRUE) THEN idx_cols_cvt  =  IMSL_LONG(idx_cols)
   IF (ARG_PRESENT(nmissing)) THEN nmissing_spc = IMSL_0
   IF (ARG_PRESENT(gcounts)) THEN gcounts_spc = IMSL_LONARR(ng_cvt)
   ; Floating point arguments and keywords
   ;
   IF (type EQ TYP_DOUBLE) THEN BEGIN
      result = dblarr(n, n)
      x_cvt = double(transpose(x))
      IF (ARG_PRESENT(sum_weights)) THEN sw_spc = dblarr(ng_cvt)
      IF (ARG_PRESENT(means)) THEN means_spc = dblarr(n, ng_cvt)
      IF (ARG_PRESENT(u)) THEN u_spc = dblarr(n, n)
   END ELSE BEGIN
      result = fltarr(n, n)
      x_cvt = float(transpose(x))
      IF (ARG_PRESENT(sum_weights)) THEN sw_spc = fltarr(ng_cvt)
      IF (ARG_PRESENT(means)) THEN means_spc = fltarr(n, ng_cvt)
      IF (ARG_PRESENT(u)) THEN u_spc = fltarr(n, n)
   END
   ;
   ; Call the system function.
   ;
   err_status = 0L
   MATHSTAT_271, type, err_status, X_cvt, x_col_dim, m, n, ng_cvt, $
                              gcounts_spc, $
                              sw_spc, $
                              nmissing_spc, $
                              means_spc, $
                              u_spc, $
                              idx_vars_cvt, $
                              idx_cols_cvt, $
                              result

   IF (ARG_PRESENT(gcounts)) THEN gcounts = gcounts_spc
   IF (ARG_PRESENT(sum_weights)) THEN sum_weights = sw_spc
   IF (ARG_PRESENT(nmissing)) THEN nmissing   =   nmissing_spc
   IF (ARG_PRESENT(means)) THEN means  =  TRANSPOSE(means_spc)
   IF (ARG_PRESENT(u)) THEN u = TRANSPOSE(u_spc)

; Return
   RETURN, TRANSPOSE(result)
END
