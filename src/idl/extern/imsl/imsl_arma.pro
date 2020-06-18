; $Id: //depot/idl/IDL_70/idl_src/libs/imsl/imsl_6.0/lib/imsl_arma.pro#1 $
;
; Copyright (c) 1970-2006, VISUAL NUMERICS Inc. All Rights Reserved.
; This software is confidential information which is proprietary to and a
; trade secret of Visual Numerics Inc.  Use, duplication or disclosure is
; subject to the terms of an appropriate license agreement.
;
FUNCTION imsl_arma, z, $                          ;INPUT 1-D array: floating point
               p, $                               ;INPUT Scalar LONG
               q, $                               ;INPUT Scalar LONG
               backward_origin=backward_origin, $ ;INPUT Scalar LONG
               itmax=itmax, $                     ;INPUT Scalar LONG
               lgth_backcast=lgth_backcast, $     ;INPUT Scalar LONG
               n_predict=n_predict, $             ;INPUT Scalar LONG
               constant=constant, $               ;INPUT Scalar ON/OFF flag
               double=double, $                   ;INPUT Scalar ON/OFF flag
               lsq=lsq, $                         ;INPUT Scalar ON/OFF flag
               moments=moments, $                 ;INPUT Scalar ON/OFF flag
               no_constant=no_constant, $         ;INPUT Scalar ON/OFF flag
               ar_lags=ar_lags, $                 ;INPUT 1-D array: LONG
               ma_lags=ma_lags, $                 ;INPUT 1-D array: LONG
               confidence=confidence, $           ;INPUT Scalar floating point
               err_rel=err_rel, $                 ;INPUT Scalar floating point
               mean_est=mean_est, $               ;INPUT Scalar floating point
               tol_backcast=tol_backcast, $       ;INPUT Scalar floating point
               tol_convergence=tol_convergence, $ ;INPUT Scalar floating point
               init_est_ar=init_est_ar, $         ;INPUT 1-D array: floating point
               init_est_ma=init_est_ma, $         ;INPUT 1-D array: floating point
               autocov=autocov, $                 ;OUTPUT 1-D array: floating point
               forecast=forecast, $               ;OUTPUT 2-D array: floating point
               param_est_cov=param_est_cov, $     ;OUTPUT 2-D array: floating point
               residual=residual, $               ;OUTPUT 1-D array: floating point
               ss_residual=ss_residual            ;OUTPUT Scalar floating point

@imsl_init.pro
   ON_ERROR, on_err_action

   ;
   ; Error checking.  
   ;   - Make sure the first input argument is a simple array. 
   ;   - Second and third arguments must be nonnegative and add up
   ;     to at least two.
   ;   - INIT_EST_AR and INIT_EST_MA must be used together.
   ;   - LGTH_BACKCAST and TOL_BACKCAST must be used together.
   ;   - INIT_EST_AR must be a 1-D array of length p
   ;   - INIT_EST_MA must be a 1-D array of length q
   ;   - AR_LAGS must be a 1-D array of length p
   ;   - MA_LAGS must be a 1-D array of length q
   ;   - If FORECAST is present, then N_PREDICT must also be present.
   ;   - If CONFIDENCE is present, then FORECAST must also be present.
   ;   - If BACKWARD_ORIGIN is present, then FORECAST must also be present.
   ;
   nargs = n_params()
   IF (nargs NE 3) THEN $
         message, "Incorrect number of arguments."
   ;
   size_z = IMSL_SIZE(z)
   IF (size_z(0) NE 1) THEN BEGIN 
      message, "Z must be a 1-D array."
   END
   nobs = size_z(1) ;Define nobs
   p_cvt = (IMSL_LONG(p))(0)
   q_cvt = (IMSL_LONG(q))(0)
   IF (((p_cvt LT 0) OR (q_cvt LT 0)) OR ((p_cvt + q_cvt) LT 1)) THEN $
     message, "The values for P and Q must be nonnegative, and their sum must be positive"
   ;
   i_tmp = KEYWORD_SET(init_est_ar) +  KEYWORD_SET(init_est_ma)
   IF ((i_tmp NE 0) AND (i_tmp NE 2)) THEN $
     message, "The keywords INIT_EST_AR and INIT_EST_MA must be used together"
   ;
   i_tmp = KEYWORD_SET(lgth_backcast) +  KEYWORD_SET(tol_backcast)
   IF ((i_tmp NE 0) AND (i_tmp NE 2)) THEN $
     message, "The keywords LGTH_BACKCAST and TOL_BACKCAST must be used together"
   ;
   IF (KEYWORD_SET(init_est_ar) EQ true) THEN BEGIN
      size_tmp = IMSL_SIZE(init_est_ar)
      IF (size_tmp(0) NE 1) THEN $
      message, "INIT_EST_AR must be a 1-D array."
      IF (N_ELEMENTS(init_est_ar) NE p_cvt) THEN $
        message, "INIT_EST_AR is not the correct length."
   END
   ;
   IF (KEYWORD_SET(init_est_ma) EQ true) THEN BEGIN
      size_tmp = IMSL_SIZE(init_est_ma)
      IF (size_tmp(0) NE 1) THEN $
      message, "INIT_EST_MA must be a 1-D array."
      IF (N_ELEMENTS(init_est_ma) NE q_cvt) THEN $
        message, "INIT_EST_MA is not the correct length."
   END
   ;
   IF (KEYWORD_SET(ar_lags) EQ true) THEN BEGIN
      size_tmp = IMSL_SIZE(ar_lags)
      IF (size_tmp(0) NE 1) THEN $
      message, "AR_LAGS must be a 1-D array."
      IF (N_ELEMENTS(ar_lags) NE p_cvt) THEN $
        message, "AR_LAGS is not the correct length."
   END
   ;
   IF (KEYWORD_SET(ma_lags) EQ true) THEN BEGIN
      size_tmp = IMSL_SIZE(ma_lags)
      IF (size_tmp(0) NE 1) THEN $
      message, "MA_LAGS must be a 1-D array."
      IF (N_ELEMENTS(ma_lags) NE q_cvt) THEN $
        message, "MA_LAGS is not the correct length."
   END
   ;
   i_tmp = ARG_PRESENT(forecast) + KEYWORD_SET(n_predict)
   IF (i_tmp EQ 1) THEN $
     message, "The keywords FORECAST and N_PREDICT must be used together."
   ;
   IF ((KEYWORD_SET(confidence) EQ TRUE) AND (ARG_PRESENT(forecast) EQ FALSE)) THEN $
     message, "The keyword CONFIDENCE can be used only if FORECAST is also present"
   ;
   IF ((KEYWORD_SET(backward_origin) EQ TRUE) AND (ARG_PRESENT(forecast) EQ FALSE)) THEN $
     message, "The keyword BACKWARD_ORIGIN can be used only if FORECAST is also present"
   ;
   ;ERROR CHECKING COMPLETE.
   ;
   ; Decide on what precision to use.
   ; Use the highest precision of the input argument(s).
   ;
   type = TYP_FLOAT
   IF (size_z(N_ELEMENTS(size_z)-2) GT TYP_FLOAT) THEN type = TYP_DOUBLE
   IF (KEYWORD_SET(double) EQ true) THEN type = TYP_DOUBLE
   ;
   ; Setup the parameters for the call to the system function.
   ;
   ; Input LONG keyword(s)
   IF (KEYWORD_SET(backward_origin) EQ TRUE) THEN $
     bkwd_origin_cvt = (IMSL_LONG(backward_origin))(0)
   IF (KEYWORD_SET(itmax) EQ TRUE) THEN $
     itmax_cvt = (IMSL_LONG(itmax))(0)
   IF (KEYWORD_SET(lgth_backcast) EQ TRUE) THEN $
     lgth_bkcst_cvt = (IMSL_LONG(lgth_backcast))(0)
   IF (KEYWORD_SET(n_predict) EQ TRUE) THEN $
     n_predict_cvt = (IMSL_LONG(n_predict))(0)
   IF (KEYWORD_SET(constant) EQ TRUE) THEN $
     constant_cvt = IMSL_1
   IF (KEYWORD_SET(lsq) EQ TRUE) THEN $
     lsq_cvt = IMSL_1
   IF (KEYWORD_SET(moments) EQ TRUE) THEN $
     moments_cvt = IMSL_1
   IF (KEYWORD_SET(no_constant) EQ TRUE) THEN $
     no_constant_cvt = IMSL_1
   IF (KEYWORD_SET(ar_lags) EQ TRUE) THEN $
     ar_lags_cvt = IMSL_LONG(ar_lags)
   IF (KEYWORD_SET(ma_lags) EQ TRUE) THEN $
     ma_lags_cvt = IMSL_LONG(ma_lags)
   ;
   ; Floating point arguments and keywords
   IF (type EQ TYP_DOUBLE) THEN BEGIN
      ; Result vector.
      result = dblarr(p_cvt+q_cvt+1)
      ; 
      ; Input
      z_cvt = double(z)
      IF (KEYWORD_SET(confidence) EQ TRUE) THEN $
        confidence_cvt = (double(confidence))(0)
      IF (KEYWORD_SET(err_rel) EQ TRUE) THEN $
        err_rel_cvt = (double(err_rel))(0)
      IF (KEYWORD_SET(mean_est) EQ TRUE) THEN $
        mean_est_cvt = (double(mean_est))(0)
      IF (KEYWORD_SET(tol_backcast) EQ TRUE) THEN $
        tol_bkcst_cvt = (double(tol_backcast))(0)
      IF (KEYWORD_SET(tol_convergence) EQ TRUE) THEN $
        tol_conv_cvt = (double(tol_convergence))(0)
      IF (KEYWORD_SET(init_est_ar) EQ TRUE) THEN $
        init_est_ar_cvt = double(init_est_ar)
      IF (KEYWORD_SET(init_est_ma) EQ TRUE) THEN $
        init_est_ma_cvt = double(init_est_ma)
      ; 
      ; Output 
      IF (ARG_PRESENT(autocov) EQ TRUE) THEN $
         autocov_spc = dblarr(p_cvt+q_cvt+2)
      IF (ARG_PRESENT(forecast) EQ TRUE) THEN BEGIN
         IF (KEYWORD_SET(backward_origin) EQ true) $
           THEN i_tmp = (bkwd_origin_cvt > IMSL_0) ELSE i_tmp = IMSL_0
         forecast_spc = dblarr(i_tmp+3, (n_predict_cvt > IMSL_0))
      END
      IF (ARG_PRESENT(param_est_cov) EQ TRUE) THEN BEGIN
         IF (KEYWORD_SET(no_constant) EQ TRUE) THEN i_tmp = p+q ELSE i_tmp = p+q+1
         param_est_spc = dblarr(i_tmp, i_tmp)
      END
      IF (ARG_PRESENT(residual) EQ TRUE) THEN BEGIN
         IF (KEYWORD_SET(ar_lags) EQ TRUE) $
           THEN max_ar_lags = IMSL_LONG(max(ar_lags_cvt)) ELSE max_ar_lags = p_cvt
         IF (KEYWORD_SET(lgth_backcast) EQ TRUE) $
           THEN length = lgth_bkcst_cvt ELSE length = 10
         residual_spc = dblarr(nobs - max_ar_lags + length)
      END
      IF (ARG_PRESENT(ss_residual) EQ TRUE) THEN $
         ss_residual_spc = double(0.0)
   END ELSE BEGIN
      ; Result vector.
      result = fltarr(p_cvt+q_cvt+1)
      ; 
      ; Input
      z_cvt = float(z)
      IF (KEYWORD_SET(confidence) EQ TRUE) THEN $
        confidence_cvt = (float(confidence))(0)
      IF (KEYWORD_SET(err_rel) EQ TRUE) THEN $
        err_rel_cvt = (float(err_rel))(0)
      IF (KEYWORD_SET(mean_est) EQ TRUE) THEN $
        mean_est_cvt = (float(mean_est))(0)
      IF (KEYWORD_SET(tol_backcast) EQ TRUE) THEN $
        tol_bkcst_cvt = (float(tol_backcast))(0)
      IF (KEYWORD_SET(tol_convergence) EQ TRUE) THEN $
        tol_conv_cvt = (float(tol_convergence))(0)
      IF (KEYWORD_SET(init_est_ar) EQ TRUE) THEN $
        init_est_ar_cvt = float(init_est_ar)
      IF (KEYWORD_SET(init_est_ma) EQ TRUE) THEN $
        init_est_ma_cvt = float(init_est_ma)
      ; 
      ; Output 
      IF (ARG_PRESENT(autocov) EQ TRUE) THEN $
         autocov_spc = fltarr(p_cvt+q_cvt+2)
      IF (ARG_PRESENT(forecast) EQ TRUE) THEN BEGIN
         IF (KEYWORD_SET(backward_origin) EQ true) $
           THEN i_tmp = (bkwd_origin_cvt > IMSL_0) ELSE i_tmp = IMSL_0
         forecast_spc = fltarr(i_tmp+3, (n_predict_cvt > IMSL_0))
      END
      IF (ARG_PRESENT(param_est_cov) EQ TRUE) THEN BEGIN
         IF (KEYWORD_SET(no_constant) EQ TRUE) THEN i_tmp = p+q ELSE i_tmp = p+q+1
         param_est_spc = fltarr(i_tmp, i_tmp)
      END
      IF (ARG_PRESENT(residual) EQ TRUE) THEN BEGIN
         IF (KEYWORD_SET(ar_lags) EQ TRUE) $
           THEN max_ar_lags = IMSL_LONG(max(ar_lags_cvt)) ELSE max_ar_lags = p_cvt
         IF (KEYWORD_SET(lgth_backcast) EQ TRUE) $
           THEN length = lgth_bkcst_cvt ELSE length = 10
         residual_spc = fltarr(nobs - max_ar_lags + length)
      END
      IF (ARG_PRESENT(ss_residual) EQ TRUE) THEN $
         ss_residual_spc = float(0.0)
   END
   ;
   ; Call the system function.
   ;
   err_status = 0L
   MATHSTAT_103,   type, err_status, $
               z_cvt, $
               p_cvt, $
               q_cvt, $
               nobs, $
               bkwd_origin_cvt, $
               itmax_cvt, $
               lgth_bkcst_cvt, $
               n_predict_cvt, $
               constant_cvt, $
               lsq_cvt, $
               moments_cvt, $
               no_constant_cvt, $
               ar_lags_cvt, $
               ma_lags_cvt, $
               confidence_cvt, $
               err_rel_cvt, $
               mean_est_cvt, $
               tol_bkcst_cvt, $
               tol_conv_cvt, $
               init_est_ar_cvt, $
               init_est_ma_cvt, $
               autocov_spc, $
               forecast_spc, $
               param_est_spc, $
               residual_spc, $
               ss_residual_spc, $
               result
                 
   ;
   ; Now copy over all output keywords results.
   ;
      IF (ARG_PRESENT(autocov) EQ TRUE) THEN $
         autocov = autocov_spc
      IF (ARG_PRESENT(forecast) EQ TRUE) THEN BEGIN
         forecast = transpose(forecast_spc)
      END
      IF (ARG_PRESENT(param_est_cov) EQ TRUE) THEN BEGIN
         param_est_cov =transpose(param_est_spc)
      END
      IF (ARG_PRESENT(residual) EQ TRUE) THEN BEGIN
         residual=residual_spc
      END
      IF (ARG_PRESENT(ss_residual) EQ TRUE) THEN $
         ss_residual=ss_residual_spc
   ;
   ; Return.
   ;
   RETURN, result
END
   

