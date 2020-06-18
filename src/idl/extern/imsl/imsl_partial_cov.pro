; $Id: //depot/idl/IDL_70/idl_src/libs/imsl/imsl_6.0/lib/imsl_partial_cov.pro#1 $
;
; Copyright (c) 1970-2006, VISUAL NUMERICS Inc. All Rights Reserved.
; This software is confidential information which is proprietary to and a
; trade secret of Visual Numerics Inc.  Use, duplication or disclosure is
; subject to the terms of an appropriate license agreement.
;
FUNCTION imsl_partial_cov, n_independent, $ ;INPUT Scalar LONG
                      n_dependent, $        ;INPUT Scalar LONG
                      X, $                  ;INPUT 2-D array: floating point 
           indices=indices, $               ;INPUT 1-D array: LONG
           cov=cov, $                       ;INPUT Scalar ON/OFF flag
           corr=corr, $                     ;INPUT Scalar ON/OFF flag
           df=df, $                         ;INPUT/OUTPUT Scalar LONG
           pvals=pvals, $                   ;OUTPUT 2-D array: floating point 
           double=double                    ;INPUT Scalar ON/OFF flag
           
@imsl_init.pro
   ON_ERROR, on_err_action
   ;
   ; Error checking:
   ;  - X must be a 2-D square array (nxn), n=n_independent+n_dependent  
   ;  - CORR and COV are mutually exclusive.
   ;  - DF and PVALS must be used together.
   ;  - If INDICES is supplied, it must be a 1D array of length n.
   ;          
   nargs = n_params()
   IF (nargs NE 3)  THEN message, 'Incorrect number of arguments.'
     
   n = IMSL_LONG(n_independent(0)+n_dependent(0))
   size_x  =  IMSL_LONG(SIZE(x))
   IF (size_x(0) NE 2) THEN $
     message, 'X must be a 2-D array.'
   IF ((size_x(1) NE size_x(2)) OR (size_x(1) NE n)) THEN $
     MESSAGE, 'X is not the correct size.'
   IF (size_x(0) EQ 2) THEN n = IMSL_LONG(size_x(2)) ELSE n = IMSL_1
   IF ((KEYWORD_SET(corr) + KEYWORD_SET(cov)) GT 1) THEN $
     message, 'The keywords CORR and COV cannot both be set.'
   IF ((ARG_PRESENT(df) + ARG_PRESENT(pvals)) EQ 1) THEN $
     MESSAGE,  'The keywords DF and PVALS must be used together.'
   ; Make sure DF is defined when both DF and PVALS are present.  Be
   ; careful since DF can be zero.
   dotest = 0
   IF ((ARG_PRESENT(df) + ARG_PRESENT(pvals)) EQ 2) THEN BEGIN
      IF (N_ELEMENTS(df) EQ 0) THEN $
        MESSAGE,    'Keyword DF must be defined.' $
        ELSE dotest = 1
   END
    
   IF (KEYWORD_SET(indices)) THEN BEGIN
      size_indices  =  IMSL_LONG(SIZE(indices))
      IF (size_indices(0) NE 1) THEN $
        message, 'INDICES must be a 1-D array.'
      IF (size_indices(1) NE n)  THEN $
        MESSAGE, 'INDICES is not the correct size.'
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
   n_dependent_cvt = IMSL_LONG(n_dependent(0))
   n_independent_cvt = IMSL_LONG(n_independent(0))
   IF (KEYWORD_SET(corr)) THEN corr_cvt = IMSL_1 ELSE cov_cvt = IMSL_1
   IF (dotest EQ 1) THEN df_cvt = IMSL_LONG(df(0))
   IF (KEYWORD_SET(indices)) THEN indices_cvt = IMSL_LONG(indices)
   ; Floating point arguments and keywords
   ;
   IF (type EQ TYP_DOUBLE) THEN BEGIN
      result = dblarr(n_dependent_cvt, n_dependent_cvt)
      IF (dotest EQ 1) THEN pvals_spc = dblarr(n_dependent_cvt, n_dependent_cvt)
   END ELSE BEGIN
      result = fltarr(n_dependent_cvt, n_dependent_cvt)
      IF (dotest EQ 1) THEN pvals_spc = fltarr(n_dependent_cvt, n_dependent_cvt)
   END
   x_cvt = imsl_cvt_arr(x, type)
   ;
   ; Call the system function.
   ;
   err_status = 0L
   MATHSTAT_267, type, err_status, x_cvt, n_independent_cvt, n_dependent_cvt, n, $
                              indices_cvt, $
                              corr_cvt, $
                              cov_cvt, $
                              df_cvt, $
                              pvals_spc, $
                              result
   IF (dotest EQ 1) THEN pvals = TRANSPOSE(pvals_spc)
   IF (dotest EQ 1) THEN df = df_cvt

; Return
   RETURN, result
END
