; $Id: //depot/idl/IDL_70/idl_src/libs/imsl/imsl_6.0/lib/imsl_allbest.pro#1 $
;
; Copyright (c) 1970-2006, VISUAL NUMERICS Inc. All Rights Reserved.
; This software is confidential information which is proprietary to and a
; trade secret of Visual Numerics Inc.  Use, duplication or disclosure is
; subject to the terms of an appropriate license agreement.
;
PRO imsl_allbest, x, $                       ;INPUT 2-D array; floating point
                  y, $                       ;INPUT 1-D array: floating point
                  Weights=weights, $         ;INPUT 1-D array; floating point
                  Frequencies=frequencies, $ ;INPUT 1-D array; floating point
                  Max_subset=max_subset, $   ;INPUT Scalar integer
                  Double=double, $           ;INPUT Scalar ON/OFF flag
                  Adj_r_squared=ars, $       ;INPUT Scalar ON/OFF flag
                  Mallows_cp=mallows_cp, $   ;INPUT Scalar ON/OFF flag
                  Max_n_best=max_n_best, $   ;INPUT Scalar integer
                  Max_n_good=max_n_good, $   ;INPUT Scalar integer
                  Cov_nobs=cov_nobs, $       ;INPUT Scalar integer
                  Cov_input=cov_input, $     ;INPUT 2-D array; floating point
                  Idx_criterions=idx_criterions, $ ;OUTPUT 1-D array: integer
                  Criterions=criterions, $   ;OUTPUT 1-D array; floating point
                  Idx_vars=idx_vars, $       ;OUTPUT 1-D array: integer
                  Indep_vars=indep_vars, $   ;OUTPUT 1-D array; integer
                  Idx_coefs=idx_coefs, $     ;OUTPUT 1-D array: integer 
                  Coefs=coefs   ;OUTPUT 2-D array; floating point    
 
@imsl_init.pro
   ON_ERROR, on_err_action
   ;
   ; Error checking.
   ; Make sure the input arguments are 1-D or 2-D arrays.
   ;
   nargs = n_params()
   IF (nargs NE 2) THEN $
     message, "Incorrect number of arguments."
   size_x = IMSL_SIZE(x)
   size_y = IMSL_SIZE(y)
   IF (size_x(0) LT 1 OR size_x(0) GT 2) THEN $
     message, "The input array, X, must be 1-D or 2-D"
   IF (size_x(0) EQ 2) THEN BEGIN
      n_rows = size_x(1)
      n_candidate = size_x(2)
   END
   IF (size_x(0) EQ 1) THEN BEGIN
      n_rows = size_x(1)
      n_candidate = IMSL_1
   END
   IF (size_y(0) NE 1) THEN $
     message, "The input array, Y, must be 1-D"
   IF (size_y(1) NE n_rows) THEN $
     message, "The input array, Y, is not the correct length"
   i = 0
   IF (KEYWORD_SET(max_subset) EQ TRUE) THEN i = i + 1   
   IF (KEYWORD_SET(ars) EQ TRUE) THEN i = i + 1
   IF (KEYWORD_SET(mallows_cp) EQ TRUE) THEN i = i + 1
   IF (i GT 1) THEN $
     message, "The following keywords are mutually exclusive: MAX_SUBSET, ADJ_R_SQUARED, and MALLOWS_CP"
   IF (KEYWORD_SET(weights) EQ TRUE) THEN BEGIN
      size_weights = IMSL_SIZE(weights)
      IF (size_weights(0) NE 1) THEN $
        message, "The weights array must be 1-D"
      IF (size_weights(1) NE n_rows) THEN $
        message, "The length of the weights is incorrect"
   END
   IF (KEYWORD_SET(frequencies) EQ TRUE) THEN BEGIN
      size_frequencies = IMSL_SIZE(frequencies)
      IF (size_frequencies(0) NE 1) THEN $
        message, "The frequencies array must be 1-D"
      IF (size_frequencies(1) NE n_rows) THEN $
        message, "The length of the frequencies is incorrect"
   END
   i = 0
   IF (ARG_PRESENT(cov_nobs) EQ TRUE) THEN i = i + 1
   IF (ARG_PRESENT(cov_input) EQ TRUE) THEN i = i + 1
   IF (i EQ 1) THEN $
     message, "COV_NOBS and COV_INPUT must be used together"
   IF (ARG_PRESENT(cov_input) EQ TRUE) THEN BEGIN
      size_cov_input = IMSL_SIZE(cov_input)
      IF (size_cov_input(0) NE 2) THEN $
        message, "The COV_INPUT array must be 2-D"
      IF (size_cov_input(1) NE (n_candidate+1)) THEN $
        message, "The COV_INPUT array is not of the correct order"
      IF (size_cov_input(2) NE (n_candidate+1)) THEN $
        message, "The COV_INPUT array is not of the correct order"
   END
   i = 0
   IF (ARG_PRESENT(idx_coefs) EQ TRUE) THEN i = i + 1
   IF (ARG_PRESENT(coefs) EQ TRUE) THEN i = i + 1
   IF (i EQ 1) THEN $
     message, "IDX_COEFS and COEFS must be used together"
   i = 0
   IF (ARG_PRESENT(idx_vars) EQ TRUE) THEN i = i + 1
   IF (ARG_PRESENT(indep_vars) EQ TRUE) THEN i = i + 1
   IF (i EQ 1) THEN $
     message, "IDX_VARS and INDEP_VARS must be used together"
   i = 0
   IF (ARG_PRESENT(idx_criterions) EQ TRUE) THEN i = i + 1
   IF (ARG_PRESENT(criterions) EQ TRUE) THEN i = i + 1
   IF (i eq 1) THEN $
      message, "IDX_CRITERIONS and CRITERIONS must be used together"
   ;
   ; Decide on what precision to use.
   ; Use the highest precision of the input argument(s).
   ;
   type = TYP_FLOAT
   IF (size_x(N_ELEMENTS(size_x)-2) GT TYP_FLOAT) THEN type = TYP_DOUBLE
   IF (size_y(N_ELEMENTS(size_y)-2) GT TYP_FLOAT) THEN type = TYP_DOUBLE
   IF (KEYWORD_SET(double) EQ true) THEN type = TYP_DOUBLE
   ; 
   ; Compute ntbest and nsize. They are used in 
   ; determining the length of some of the output keywords.
   ;
   IF (KEYWORD_SET(max_n_best) EQ TRUE) THEN BEGIN
      max_n_best_local = (IMSL_LONG(max_n_best))(0)
   END ELSE BEGIN
      max_n_best_local = IMSL_1
   END
   IF (KEYWORD_SET(max_subset) EQ TRUE) THEN BEGIN
      max_subset_local = (IMSL_LONG(max_subset))(0)
      nsize = max_subset_local
      ntbest = max_n_best_local*max_subset_local
   END ELSE IF (KEYWORD_SET(ars) EQ TRUE) THEN BEGIN
      nsize = n_candidate
      ntbest = max_n_best_local
      max_subset_local = n_candidate
   END ELSE IF (KEYWORD_SET(mallows_cp) EQ TRUE) THEN BEGIN
      nsize = n_candidate 
      ntbest = max_n_best_local
      max_subset_local = n_candidate
   END ELSE BEGIN
      nsize = n_candidate
      ntbest = max_n_best_local*n_candidate
      max_subset_local = n_candidate
   END
   ;
   ; Setup the parameters for the call to the system function.
   ;
   IF (KEYWORD_SET(max_subset) EQ TRUE) THEN BEGIN
      max_subset_cvt = (IMSL_LONG(max_subset))(0)
      IF (max_subset_cvt LT 1) THEN $
        message, "MAX_SUBSET must be positive" 
   END   
   IF (KEYWORD_SET(ars) EQ TRUE) THEN $
        ars = IMSL_1
   IF (KEYWORD_SET(mallows_cp) EQ TRUE) THEN $
        mallows_cp = IMSL_1
   IF (KEYWORD_SET(max_n_best) EQ TRUE) THEN BEGIN
      max_n_best_cvt = (IMSL_LONG(max_n_best))(0)
      IF (max_n_best_cvt LT 1) THEN $
        message, "MAX_N_BEST must be positive" 
   END
   IF (KEYWORD_SET(max_n_good) EQ TRUE) THEN BEGIN
      max_n_good_cvt = (IMSL_LONG(max_n_good))(0)
      max_n_good_local = max_n_good_cvt
      IF (max_n_good_cvt LT 1) THEN $
        message, "MAX_N_BEST must be positive"
   END ELSE BEGIN
      max_n_good_local = IMSL_LONG(10)
   END
   IF (type EQ TYP_DOUBLE) THEN BEGIN
      x_cvt = double(x)
      y_cvt = double(y)
      IF (KEYWORD_SET(weights) EQ TRUE) THEN $
        weights_cvt = double(weights)
      IF (KEYWORD_SET(frequencies) EQ TRUE) THEN $
        frequencies_cvt = double(frequencies)
      IF (KEYWORD_SET(cov_nobs) EQ TRUE) THEN $
        cov_nobs_cvt = (IMSL_LONG(cov_nobs))(0)
      IF (KEYWORD_SET(cov_input) EQ TRUE) THEN $
        cov_input_cvt = double(cov_input)
      IF (ARG_PRESENT(criterions) EQ TRUE) THEN BEGIN
         idx_criterions_spc = IMSL_LONARR(nsize+1)
         IF (max_n_good_local*nsize GT n_candidate) THEN BEGIN
            junk = max_n_good_local*nsize
         END ELSE BEGIN
            junk = n_candidate
         END
         criterions_spc = dblarr(junk)
      END
      IF (ARG_PRESENT(indep_vars) EQ TRUE) THEN BEGIN
         idx_vars_spc = IMSL_LONARR(nsize+1)
         junk = max_n_good_local*nsize*(nsize+1)/2
         indep_vars_spc = IMSL_LONARR(junk)
      END
      IF (ARG_PRESENT(coefs) EQ TRUE) THEN BEGIN
         idx_coefs_spc = IMSL_LONARR(ntbest+1)
         junk = (max_n_best_local*max_subset_local*(max_subset_local+1))/2
         coefs_spc = dblarr(5, junk)
      END
   END ELSE BEGIN
      x_cvt = float(x)
      y_cvt = float(y)
      IF (KEYWORD_SET(weights) EQ TRUE) THEN $
        weights_cvt = float(weights)
      IF (KEYWORD_SET(frequencies) EQ TRUE) THEN $
        frequencies_cvt = float(frequencies)
      IF (KEYWORD_SET(cov_nobs) EQ TRUE) THEN $
        cov_nobs_cvt = (IMSL_LONG(cov_nobs))(0)
      IF (KEYWORD_SET(cov_input) EQ TRUE) THEN $
        cov_input_cvt = float(cov_input)
      IF (ARG_PRESENT(criterions) EQ TRUE) THEN BEGIN
         idx_criterions_spc = IMSL_LONARR(nsize+1)
         IF (max_n_good_local*nsize GT n_candidate) THEN BEGIN
            junk = max_n_good_local*nsize
         END ELSE BEGIN
            junk = n_candidate
         END
         criterions_spc = fltarr(junk)
      END
      IF (ARG_PRESENT(indep_vars) EQ TRUE) THEN BEGIN
         idx_vars_spc = IMSL_LONARR(nsize+1)
         junk = max_n_good_local*nsize*(nsize+1)/2
         indep_vars_spc = IMSL_LONARR(junk)
      END
      IF (ARG_PRESENT(coefs) EQ TRUE) THEN BEGIN
         idx_coefs_spc = IMSL_LONARR(ntbest+1)
         junk = (max_n_best_local*max_subset_local*(max_subset_local+1))/2
         coefs_spc = fltarr(5, junk)
      END
   END
   ;
   ; Transpose the input matrix, X.
   ;
   x_cvt = transpose(x_cvt)
   ;
   ; Call the system function.
   ;
   err_status = 0L
   MATHSTAT_100, type, err_status, x_cvt, y_cvt, n_rows, n_candidate, $
     weights_cvt, frequencies_cvt, max_subset_cvt, ars, $
     mallows_cp, max_n_best_cvt, max_n_good_cvt, cov_nobs_cvt, $
     cov_input_cvt, idx_criterions_spc, criterions_spc, $
     idx_vars_spc, indep_vars_spc, idx_coefs_spc, coefs_spc
   
   ;
   ; Now copy over all output keywords results.
   ; Since we sent an upperbound of space needed by C/Stat for some of
   ; the output keywords, we copy them into variables of the correct 
   ; size before we return them to the user.
   ;
   IF (ARG_PRESENT(criterions) EQ TRUE) THEN BEGIN
      idx_criterions = idx_criterions_spc
      junk = idx_criterions(nsize)
      IF (n_candidate GT junk) THEN junk = n_candidate
      criterions = criterions_spc(0:junk-1)
   END
   IF (ARG_PRESENT(coefs) EQ TRUE) THEN BEGIN
      idx_coefs = idx_coefs_spc
      junk = idx_coefs(ntbest)
      coefs = coefs_spc(*, 0:junk-1)
      coefs = transpose(coefs)
   END
   IF (ARG_PRESENT(indep_vars) EQ TRUE) THEN BEGIN
      idx_vars = idx_vars_spc
      junk = idx_vars(nsize)
      indep_vars = indep_vars_spc(0:junk-1)
   END
   ;
   ; Return.
   ;
   RETURN
END
   

