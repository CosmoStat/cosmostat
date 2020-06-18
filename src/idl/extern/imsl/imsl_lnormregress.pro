; $Id: //depot/idl/IDL_70/idl_src/libs/imsl/imsl_6.0/lib/imsl_lnormregress.pro#1 $
;
; Copyright (c) 1970-2006, VISUAL NUMERICS Inc. All Rights Reserved.
; This software is confidential information which is proprietary to and a
; trade secret of Visual Numerics Inc.  Use, duplication or disclosure is
; subject to the terms of an appropriate license agreement.
;
FUNCTION imsl_lnormregress, x, $                       ;INPUT 2-D array: floating point
                   y, $                       ;INPUT 1-D array: floating point
                   lav=lav, $                 ;INPUT Scalar ON/OFF flag
                   llp=llp, $                 ;INPUT Scalar ON/OFF flag
                   lmv=lmv, $                 ;INPUT Scalar ON/OFF flag
                   p=p, $                     ;INPUT Scalar floating point
                   eps=eps, $                 ;INPUT Scalar floating point
                   weights=weights, $         ;INPUT 1-D array: floating point
                   frequencies=frequencies, $ ;INPUT 1-D array: floating point
                   no_intercept=no_intercept, $   ;INPUT Scalar ON/OFF flag
                   nmissing=nmissing, $       ;OUTPUT Scalar LONG
                   rank=rank, $               ;OUTPUT Scalar LONG
                   iters=iters, $             ;OUTPUT Scalar LONG
                   sea=sea, $                 ;OUTPUT Scalar floating point
                   scale=scale, $             ;OUTPUT Scalar floating point
                   resid_max=resid_max, $     ;OUTPUT Scalar floating point
                   df=df, $                   ;OUTPUT Scalar LONG
                   tolerance=tolerance, $     ;INPUT Scalar floating point 
                   r_matrix=r_matrix, $       ;OUTPUT 2-D array: floating point 
                   residuals=residuals, $     ;OUTPUT 1-D array: floating point 
                   resid_norm=resid_norm, $   ;OUTPUT Scalar LONG                 
                   double=DOUBLE            ;INPUT Scalar ON/OFF flag

@imsl_init.pro 
   ON_ERROR, on_err_action

   ;
   ; Error checking:
   ; The following checks are performed.
   ;   X must be a 1D or 2D array.
   ;   Y must be a 1D array.
   ;   If supplied, WEIGHTS must be a 1D array.
   ;   If supplied, FREQUENCIES must be a 1D array.
   ;   Check that certain keywords are not supplied for some methods.
   nargs = n_params()
   IF (nargs NE 2) THEN message, "Incorrect number of arguments."
   size_x = IMSL_SIZE(x)
   IF ((size_x(0) NE 1) AND (size_x(0) NE 2)) THEN BEGIN
      message, "X must be a 1-D or 2-D array."
   END
   IF (size_x(0) EQ 2) THEN BEGIN
      n_rows = size_x(1)
      n_indep = size_x(2)
   END
   IF (size_x(0) EQ 1) THEN BEGIN
      n_rows = size_x(1)
      n_indep = IMSL_1
   END
   size_y = IMSL_SIZE(y)
   IF (size_y(0) NE 1) THEN BEGIN
      message, "Y must be a 1-D array."
   END
   IF (size_y(1) NE n_rows) THEN BEGIN
      message, "Y is not the correct size."
   END
   IF (KEYWORD_SET(weights)) THEN BEGIN 
      size_weights = IMSL_SIZE(weights)
      IF (size_weights(0) NE 1) THEN message, 'WEIGHTS must be a 1-D array.'
      IF (size_weights(1) NE n_rows) THEN $
        message, 'WEIGHTS is not the correct size.'
   END
   IF (KEYWORD_SET(frequencies)) THEN BEGIN 
      size_frequencies = IMSL_SIZE(frequencies)
      IF (size_frequencies(0) NE 1) THEN message, 'FREQUENCIES must be a 1-D array.'
      IF (size_frequencies(1) NE n_rows) THEN $
        message, 'FREQUENCIES is not the correct size.'
   END
   IF ((KEYWORD_SET(lav) + KEYWORD_SET(llp) + KEYWORD_SET(lmv)) GT 1) THEN $
     MESSAGE, 'At most one method may be specified. '
   method = IMSL_1 ; default
   IF (KEYWORD_SET(llp)) THEN method = IMSL_2
   IF (KEYWORD_SET(lmv)) THEN method = IMSL_3
   IF (method NE 2) THEN BEGIN
      IF (KEYWORD_SET(weights)) THEN MESSAGE, "Keyword WEIGHTS requires the keyword LLP."
      IF (KEYWORD_SET(frequencies)) THEN MESSAGE, "Keyword FREQUENCIES requires the keyword LLP."
      IF (KEYWORD_SET(tolerance)) THEN MESSAGE, "Keyword TOLERANCE requires the keyword LLP."
      IF (ARG_PRESENT(r_matrix)) THEN MESSAGE, "Keyword R_MATRIX requires the keyword LLP."
      IF (ARG_PRESENT(eps)) THEN MESSAGE, "Keyword EPS requires the keyword LLP."
      IF (ARG_PRESENT(df)) THEN MESSAGE, "Keyword DF requires the keyword LLP."
      IF (ARG_PRESENT(residuals)) THEN MESSAGE, "Keyword RESIDUALS requires the keyword LLP."
      IF (ARG_PRESENT(scale)) THEN MESSAGE, "Keyword SCALE requires the keyword LLP."
      IF (ARG_PRESENT(resid_norm)) THEN MESSAGE, "Keyword RESID_NORM requires the keyword LLP."
   END ELSE BEGIN
      IF (NOT KEYWORD_SET(p)) THEN MESSAGE,  "Keyword LLP requires the keyword P."
   END
   IF (method NE 1) THEN BEGIN
      IF (ARG_PRESENT(sea)) THEN MESSAGE, "Keyword SEA requires the LAV method."
   END 
   IF (method NE 3) THEN BEGIN
      IF (ARG_PRESENT(resid_max)) THEN MESSAGE, "Keyword RESID_MAX requires the LMV method."
   END           
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
   IF (KEYWORD_SET(no_intercept)) THEN no_intercept_cvt = IMSL_1
   ; Output LONG keyword(s)
   IF (ARG_PRESENT(rank) EQ TRUE) THEN rank_spc = IMSL_0
   IF (ARG_PRESENT(iters) EQ TRUE) THEN iters_spc = IMSL_0
   IF (ARG_PRESENT(nmissing) EQ TRUE) THEN nmissing_spc = IMSL_0
   ;
   ; Floating point arguments and keywords
   IF (KEYWORD_SET(no_intercept)) THEN ncoef = n_indep ELSE ncoef = n_indep + 1
   IF (type EQ TYP_DOUBLE) THEN BEGIN
      mach = imsl_machine(/DOUBLE)
      ; Input
      x_cvt  =  DOUBLE(TRANSPOSE(x))
      y_cvt  =  DOUBLE(y)      
      IF (KEYWORD_SET(eps) EQ TRUE) THEN eps_cvt = DOUBLE(eps(0))
      IF (KEYWORD_SET(p) EQ TRUE) THEN p_cvt = DOUBLE(p(0))
      IF (KEYWORD_SET(eps)) THEN eps_cvt  =  DOUBLE(eps(0)) $
        ELSE eps_CVT = 100.*(mach.max_rel_space)
      IF (KEYWORD_SET(tolerance)) THEN tolerance_cvt  =  DOUBLE(tolerance(0)) $
        ELSE tolerance_CVT = 100.*(mach.max_rel_space)
      IF (KEYWORD_SET(weights) EQ TRUE) THEN weights_cvt = DOUBLE(weights)
      IF (KEYWORD_SET(frequencies) EQ TRUE) THEN frequencies_cvt = DOUBLE(frequencies)
      ; Output
      IF (ARG_PRESENT(df) EQ TRUE) THEN df_spc = DOUBLE(0.0)
      IF (ARG_PRESENT(sea) EQ TRUE) THEN sea_spc = DOUBLE(0.0)
      IF (ARG_PRESENT(resid_max) EQ TRUE) THEN resid_max_spc = DOUBLE(0.0)
      IF (ARG_PRESENT(scale) EQ TRUE) THEN scale_spc = DOUBLE(0.0)
      IF (ARG_PRESENT(resid_norm) EQ TRUE) THEN resid_norm_spc = DOUBLE(0.0)
      IF (ARG_PRESENT(residuals) EQ TRUE) THEN residuals_spc = DBLARR(n_rows)
      IF (ARG_PRESENT(r_matrix) EQ TRUE) THEN  r_matrix_spc = DBLARR(ncoef, ncoef)
      result = DBLARR(ncoef)
   END ELSE BEGIN
      mach  =  imsl_machine(/FLOAT)
      ; Input
      x_cvt  =  FLOAT(TRANSPOSE(x))
      y_cvt  =  FLOAT(y)      
      IF (KEYWORD_SET(p) EQ TRUE) THEN p_cvt = FLOAT(p(0))
      IF (KEYWORD_SET(eps)) THEN eps_cvt  =  FLOAT(eps(0)) $
        ELSE eps_CVT = 100.*(mach.max_rel_space)
      IF (KEYWORD_SET(tolerance)) THEN tolerance_cvt  =  FLOAT(tolerance(0)) $
        ELSE tolerance_CVT = 100.*(mach.max_rel_space)
      IF (KEYWORD_SET(weights) EQ TRUE) THEN weights_cvt = FLOAT(weights)
      IF (KEYWORD_SET(frequencies) EQ TRUE) THEN frequencies_cvt = FLOAT(frequencies)
      ; Output
      IF (ARG_PRESENT(df) EQ TRUE) THEN df_spc = FLOAT(0.0)
      IF (ARG_PRESENT(sea) EQ TRUE) THEN sea_spc = FLOAT(0.0)
      IF (ARG_PRESENT(resid_max) EQ TRUE) THEN resid_max_spc = FLOAT(0.0)
      IF (ARG_PRESENT(scale) EQ TRUE) THEN scale_spc = FLOAT(0.0)
      IF (ARG_PRESENT(resid_norm) EQ TRUE) THEN resid_norm_spc = FLOAT(0.0)
      IF (ARG_PRESENT(residuals) EQ TRUE) THEN residuals_spc = FLTARR(n_rows)
      IF (ARG_PRESENT(r_matrix) EQ TRUE) THEN  r_matrix_spc = FLTARR(ncoef, ncoef)
      result = FLTARR(ncoef)
   END
   ;
   ; Call the system function.
   ;
   err_status = 0L

   MATHSTAT_279,type,err_status, $
                           n_rows, $
                           n_indep, $
                           method, $
                           x_cvt, $
                           y_cvt, $
                           no_intercept_cvt, $
                           eps_cvt, $
                           p_cvt, $
                           tolerance_cvt, $
                           weights_cvt, $
                           frequencies_cvt, $
                           sea_spc, $
                           resid_max_spc, $
                           resid_norm_spc, $
                           scale_spc, $
                           residuals_spc, $
                           r_matrix_spc, $
                           rank_spc, $
                           iters_spc, $
                           df_spc, $
                           nmissing_spc, $
                           result
   ;
   ; Now copy over all output keywords results.
   ;
   IF (ARG_PRESENT(sea) EQ TRUE) THEN sea = sea_spc
   IF (ARG_PRESENT(resid_max) EQ TRUE) THEN resid_max = resid_max_spc
   IF (ARG_PRESENT(scale) EQ TRUE) THEN scale = scale_spc
   IF (ARG_PRESENT(resid_norm) EQ TRUE) THEN resid_norm = resid_norm_spc
   IF (ARG_PRESENT(residuals) EQ TRUE) THEN residuals = residuals_spc
   IF (ARG_PRESENT(r_matrix) EQ TRUE) THEN  r_matrix = TRANSPOSE(r_matrix_spc)
   IF (ARG_PRESENT(rank) EQ TRUE) THEN rank = rank_spc
   IF (ARG_PRESENT(iters) EQ TRUE) THEN iters = iters_spc
   IF (ARG_PRESENT(df) EQ TRUE) THEN df = df_spc
   IF (ARG_PRESENT(nmissing) EQ TRUE) THEN nmissing = nmissing_spc
   ;
   ; Return.
   ;
   RETURN, result
END

                   
                   
                   

  
      

  
