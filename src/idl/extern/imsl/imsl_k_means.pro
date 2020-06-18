; $Id: //depot/idl/IDL_70/idl_src/libs/imsl/imsl_6.0/lib/imsl_k_means.pro#1 $
;
; Copyright (c) 1970-2006, VISUAL NUMERICS Inc. All Rights Reserved.
; This software is confidential information which is proprietary to and a
; trade secret of Visual Numerics Inc.  Use, duplication or disclosure is
; subject to the terms of an appropriate license agreement.
;
FUNCTION imsl_k_means, x, $                        ;INPUT 1-D or 2-D array: floating point
                  seeds, $                         ;INPUT 1-D or 2-D array: floating point
                  Var_Columns=var_columns, $       ;INPUT 1-D array: LONG
                  Weights=weights, $               ;INPUT 1-D array: floating point
                  Frequencies=frequencies, $       ;INPUT 1-D array: floating point
                  Itmax=itmax, $                   ;INPUT Scalar LONG
                  Double=double, $                 ;INPUT Scalar ON/OFF flag
                  Means_Cluster=means_cluster, $   ;OUTPUT 2-D array: floating point
                  Ssq_Cluster=ssq_cluster, $       ;OUTPUT 1-D array: floating point
                  Counts_Cluster=counts_cluster    ;OUTPUT 1-D array: floating point
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
   IF ((size_x(0) LT 1) OR (size_x(0) GT 2)) THEN BEGIN 
      message, "X, the array containing the observations to be clustered, must be either 1-D or 2-D"
   END
   nobs = size_x(1)
   IF (size_x(0) EQ 1) THEN BEGIN
      nvar = IMSL_1
   END ELSE BEGIN
      nvar = size_x(2)
   END
   IF (KEYWORD_SET(weights) EQ TRUE) THEN BEGIN
     size_weights = IMSL_SIZE(weights)
     IF (size_weights(0) NE 1) THEN BEGIN
        message, "The WEIGHTS array must be 1-D"
     END
     IF (size_weights(N_ELEMENTS(size_weights)-1) NE nobs) THEN BEGIN
        message, "The length of the WEIGHTS array be equal to the first dimension of X"
     END
   END
   IF (KEYWORD_SET(frequencies) EQ TRUE) THEN BEGIN
     size_frequencies = IMSL_SIZE(frequencies)
     IF (size_frequencies(0) NE 1) THEN BEGIN
        message, "The FREQUENCIES array must be 1-D"
     END
     IF (size_frequencies(N_ELEMENTS(size_frequencies)-1) NE nobs) THEN BEGIN
        message, "The length of the FREQUENCIES array must be equal to the first dimension of X"
     END
   END
   x_col_dim = IMSL_LONG(nvar)
   IF (KEYWORD_SET(var_columns) EQ TRUE) THEN BEGIN
     size_tmp = IMSL_SIZE(var_columns)
     IF (size_tmp(0) NE 1) THEN BEGIN
        message, "The VAR_COLUMNS array must be 1-D"
     END
     nvar = size_tmp(N_ELEMENTS(size_tmp)-1)
     IF (nvar GT x_col_dim) THEN BEGIN
        message, "The length of the VAR_COLUMNS array must be less than " + $
                 "or equal to the second dimension of X"
     END
   END
   size_seeds = IMSL_SIZE(seeds)
   IF ((size_seeds(0) LT 1) OR (size_seeds(0) GT 2)) THEN BEGIN
      message, "SEEDS, the array containing the cluster seeds, must be either 1-D or 2-D"
   END
   seeds_dim0 = size_seeds(1)
   n_clusters = seeds_dim0
   IF (size_seeds(0) EQ 1) THEN BEGIN
      seeds_dim1 = IMSL_1
   END ELSE BEGIN
      seeds_dim1 = size_seeds(2)
   END
   IF (seeds_dim1 NE nvar) THEN BEGIN
      IF (KEYWORD_SET(var_columns) EQ TRUE) THEN  $
        message, "The second dimension of SEEDS must be equal to the length of VAR_COLUMNS." $
        ELSE message, "The second dimension of SEEDS must be equal to the second dimension of X."
   END
   ;
   ; Decide on what precision to use.
   ; Use the highest precision of the input argument(s).
   ;
   type = TYP_FLOAT
   IF (size_x(N_ELEMENTS(size_x)-2) GT TYP_FLOAT) THEN type = TYP_DOUBLE
   IF (size_seeds(N_ELEMENTS(size_seeds)-2) GT TYP_FLOAT) THEN type = TYP_DOUBLE
   IF (KEYWORD_SET(double) EQ true) THEN type = TYP_DOUBLE
   ;
   ; Setup the parameters for the call to the system function.
   ;
   ; Result array
   l_result = IMSL_LONARR(nobs)
   ;
   ; Input LONG Keywords
   IF (KEYWORD_SET(var_columns) EQ TRUE) THEN $
        var_columns_cvt = IMSL_LONG(var_columns)
   IF (KEYWORD_SET(itmax) EQ TRUE) THEN $
     itmax_cvt = (IMSL_LONG(itmax))(0)
   ;
   ; Output LONG Keywords
   IF (ARG_PRESENT(counts_cluster) EQ TRUE) THEN $
     counts_cl_spc = IMSL_LONARR(n_clusters)
   ;
   ; Floating point arguments and keywords
   IF (type EQ TYP_DOUBLE) THEN BEGIN
      ; 
      ; Input 
      x_cvt = double(transpose(x))
      seeds_cvt = double(transpose(seeds))
      IF (KEYWORD_SET(weights) EQ TRUE) THEN $
        weights_cvt = double(weights)
      IF (KEYWORD_SET(frequencies) EQ TRUE) THEN $
        frequencies_cvt = double(frequencies)
      ; 
      ; Output 
      IF (ARG_PRESENT(means_cluster) EQ TRUE) THEN $
        means_clstr_spc = dblarr(nvar, n_clusters)
      IF (ARG_PRESENT(ssq_cluster) EQ TRUE) THEN $
        ssq_cluster_spc = dblarr(n_clusters)
   END ELSE BEGIN
      ; 
      ; Input 
      x_cvt = float(transpose(x))
      seeds_cvt = float(transpose(seeds))
      IF (KEYWORD_SET(weights) EQ TRUE) THEN $
        weights_cvt = float(weights)
      IF (KEYWORD_SET(frequencies) EQ TRUE) THEN $
        frequencies_cvt = float(frequencies)
      ; 
      ; Output 
      IF (ARG_PRESENT(means_cluster) EQ TRUE) THEN $
        means_clstr_spc = fltarr(nvar, n_clusters)
      IF (ARG_PRESENT(ssq_cluster) EQ TRUE) THEN $
        ssq_cluster_spc = fltarr(n_clusters)
   END
   ;
   ; Call the system function.
   ;
    err_status = 0L
  MATHSTAT_153, type,  err_status, x_cvt, x_col_dim, nobs, nvar, seeds_cvt, n_clusters, $
                weights_cvt, $
                frequencies_cvt, $
                var_columns_cvt, $
                itmax_cvt, $
                means_clstr_spc, $
                ssq_cluster_spc, $
                counts_cl_spc, $
                l_result
   ;
   ; Now copy over all output keywords results.
   ;
   IF (ARG_PRESENT(counts_cluster) EQ TRUE) THEN $
     counts_cluster=counts_cl_spc
      IF (ARG_PRESENT(means_cluster) EQ TRUE) THEN $
        means_cluster=transpose(means_clstr_spc) ; NOTE: transpose()
      IF (ARG_PRESENT(ssq_cluster) EQ TRUE) THEN $
        ssq_cluster=ssq_cluster_spc
   ;
   ; Return.
   ;
   RETURN, l_result
END
   

