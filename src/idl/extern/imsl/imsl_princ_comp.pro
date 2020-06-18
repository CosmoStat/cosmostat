; $Id: //depot/idl/IDL_70/idl_src/libs/imsl/imsl_6.0/lib/imsl_princ_comp.pro#1 $
;
; Copyright (c) 1970-2006, VISUAL NUMERICS Inc. All Rights Reserved.
; This software is confidential information which is proprietary to and a
; trade secret of Visual Numerics Inc.  Use, duplication or disclosure is
; subject to the terms of an appropriate license agreement.
;
FUNCTION imsl_princ_comp, covariances, $          ;INPUT 2-D array: floating point
                     Cov_Matrix=cov_matrix, $     ;INPUT Scalar ON/OFF flag
                     Corr_Matrix=corr_matrix, $   ;INPUT Scalar ON/OFF flag
                     Df=df, $                     ;INPUT Scalar LONG
                     Double=double, $             ;INPUT Scalar ON/OFF fl
                     Stdev=stdev, $               ;OUTPUT 1-D array: floating point
                     Cum_Percent=cum_percent, $   ;OUTPUT 1-D array: floating point
                     Eigenvectors=eigenvectors, $ ;OUTPUT 2-D array: floating point
                     Correlations=correlations    ;OUTPUT 2-D array: floating point
 
@imsl_init.pro
   ON_ERROR, on_err_action
   ;
   ; Error checking.
   ; Make sure the input arguments are 1-D or 2-D arrays.
   ;
   nargs = n_params()
   IF (nargs NE 1) THEN $
         message, "Incorrect number of arguments."
   size_c = IMSL_SIZE(covariances)
   IF (size_c(0) NE 2) THEN BEGIN 
      message, "The input array, COVARIANCES, must be 2-D"
   END
   IF (size_c(1) NE size_c(2)) THEN BEGIN 
      message, "The input array, COVARIANCES, must be square"
   END
   nvar = size_c(1)
   IF ((KEYWORD_SET(df) AND (NOT ARG_PRESENT(stdev))) OR $
       ((NOT KEYWORD_SET(df)) AND ARG_PRESENT(stdev)) $
      ) THEN BEGIN
        message, "The keywords STDEV and DF must be used together"
   END
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
   ; Input LONG keyword(s)
   IF (KEYWORD_SET(df) EQ TRUE) THEN $
     df_cvt = IMSL_LONG(df)
   ;
   ; Floating point arguments and keywords
   IF (type EQ TYP_DOUBLE) THEN BEGIN
      ; Result vector.
      l_result = dblarr(nvar)
      ; 
      ; Input 
      covariances_cvt = double(covariances)
      IF (KEYWORD_SET(cov_matrix) EQ TRUE) THEN $
        cov_matrix_cvt = IMSL_1
      IF (KEYWORD_SET(corr_matrix) EQ TRUE) THEN $
        corr_matrix_cvt = IMSL_1
      ; 
      ; Output 
      IF (ARG_PRESENT(cum_percent) EQ TRUE) THEN $
        cum_percent_spc = dblarr(nvar)
      IF (ARG_PRESENT(eigenvectors) EQ TRUE) THEN $
        eigenvectrs_spc = dblarr(nvar, nvar)
      IF (ARG_PRESENT(correlations) EQ TRUE) THEN $
        correlats_spc = dblarr(nvar, nvar)
      IF (ARG_PRESENT(stdev) EQ TRUE) THEN $
        stdev_spc = dblarr(nvar)
   END ELSE BEGIN
      ; Result vector.
      l_result = fltarr(nvar)
      ; 
      ; Input 
      covariances_cvt = float(covariances)
      IF (KEYWORD_SET(cov_matrix) EQ TRUE) THEN $
        cov_matrix_cvt = IMSL_1
      IF (KEYWORD_SET(corr_matrix) EQ TRUE) THEN $
        corr_matrix_cvt = IMSL_1
      ; 
      ; Output 
      IF (ARG_PRESENT(cum_percent) EQ TRUE) THEN $
        cum_percent_spc = fltarr(nvar)
      IF (ARG_PRESENT(eigenvectors) EQ TRUE) THEN $
        eigenvectrs_spc = fltarr(nvar, nvar)
      IF (ARG_PRESENT(correlations) EQ TRUE) THEN $
        correlats_spc = fltarr(nvar, nvar)
      IF (ARG_PRESENT(stdev) EQ TRUE) THEN $
        stdev_spc = fltarr(nvar)
   END
   ;
   ; Call the system function.
   ;
   err_status = 0L
   MATHSTAT_175, type,  err_status,covariances_cvt, nvar, $
     cov_matrix_cvt, corr_matrix_cvt,  cum_percent_spc, $
     eigenvectrs_spc, correlats_spc, df_cvt, stdev_spc, l_result
   ;
   ; Now copy over all output keywords results.
   ;
   IF (ARG_PRESENT(correlations) EQ TRUE) THEN $
     correlations = transpose(correlats_spc) ;NOTE transpose()
   IF (ARG_PRESENT(cum_percent) EQ TRUE) THEN $
     cum_percent = cum_percent_spc
   IF (ARG_PRESENT(eigenvectors) EQ TRUE) THEN $
     eigenvectors = transpose(eigenvectrs_spc) ;NOTE transpose()
   IF (ARG_PRESENT(stdev) EQ TRUE) THEN $
     stdev = stdev_spc
   ;
   ; Return.
   ;
   RETURN, l_result
END
   

