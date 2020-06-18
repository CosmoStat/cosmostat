; $Id: //depot/idl/IDL_70/idl_src/libs/imsl/imsl_6.0/lib/imsl_survival_glm.pro#1 $
;
; Copyright (c) 1970-2006, VISUAL NUMERICS Inc. All Rights Reserved.
; This software is confidential information which is proprietary to and a
; trade secret of Visual Numerics Inc.  Use, duplication or disclosure is
; subject to the terms of an appropriate license agreement.
;
FUNCTION imsl_survival_glm, $  
                     n_class, $                  ;INPUT Scalar LONG
                     n_continuous, $             ;INPUT Scalar LONG
                     model, $                    ;INPUT Scalar LONG
                     x, $                        ;INPUT 2-D array: LONG
                     icen=icen, $                ;INPUT Scalar LONG
                     irt=irt, $                  ;INPUT Scalar LONG
                     ilt=ilt, $                  ;INPUT Scalar LONG
                     ifreq=ifreq, $              ;INPUT Scalar LONG
                     ifix=ifix, $                ;INPUT Scalar LONG
                     eps=eps, $                  ;INPUT Scalar floating point
                     itmax=itmax, $              ;INPUT Scalar LONG
                     no_intercept=no_intercept, $;INPUT Scalar LONG
                     lp_max=lp_max, $            ;INPUT Scalar LONG
                     indices_effects=indices_effects, $;INPUT 1-D array: LONG
                     var_effects=var_effects, $  ;INPUT 1-D array: LONG
                     init_est=init_est, $        ;INPUT 1-D array: floating point
                     max_class=max_class, $      ;INPUT Scalar LONG
                     n_class_vals=n_class_vals, $;OUTPUT 1-D array: LONG
                     class_vals=class_vals, $    ;OUTPUT 1-D array: floating point
                     coef_stat=coef_stat, $      ;OUTPUT 2-D array: floating point
                     criterion=criterion, $      ;OUTPUT Scalar floating point
                     covariances=covariances, $  ;OUTPUT 2-D array: floating point
                     means=means, $              ;OUTPUT 1-D array: floating point
                     case_analysis=case_analysis, $ ;OUTPUT 2-D array: floating point
                     last_step=last_step, $      ;OUTPUT 1-D array: floating point
                     obs_status=obs_status, $    ;OUTPUT 1-D array: floating point
                     iterations=iterations, $    ;OUTPUT 2-D array: floating point
                     double=double, $            ;INPUT Scalar ON/OFF flag
                     est_nobs=est_nobs, $        ;INPUT Scalar LONG
                     est_time=est_time, $        ;INPUT Scalar LONG
                     est_npt=est_npt, $          ;INPUT Scalar LONG
                     est_delta=est_delta, $      ;INPUT Scalar LONG
                     est_prob=est_prob, $        ;OUTPUT 2-D array: floating point
                     est_xbeta=est_xbeta, $      ;OUTPUT 1-D array: floating point
                     nmissing=nmissing           ;OUTPUT Scalar LONG

@imsl_init.pro
   ON_ERROR, on_err_action

   ;
   ; Error checking:
   ; The following checks are performed.
   ;   X must be a 2D array (nobs by n_var+n_continuous + m)
   ;   VAR_EFFECTS and INDICES_EFFECTS must be used together, and both
   ;      must be 1D arrays of lengths (n_effects &  total(var_effects))
   ;   If supplied, INIT_EST must be a 1D array. (n_coef_input)
   ;   Keywords est_* must be used together.
   ;
   nargs = n_params()
   IF (nargs NE 4) THEN message, "Incorrect number of arguments."
   size_x = IMSL_SIZE(x)
   IF (size_x(0) NE 2) THEN BEGIN
      message, "X must be a 2-D array."
   END
   nobs_cvt   =   size_x(1)
   x_col_dim   =   size_x(2)
   n_class_cvt  =  IMSL_LONG(n_class(0))
   n_continuous_cvt  =  IMSL_LONG(n_continuous(0))
   m_cvt  =  size_x(2) - n_class_cvt - n_continuous_cvt
   IF (m_cvt LT 1) THEN $
     MESSAGE,   "The second dimension of X is not large enough."
   ; Be careful in case IRT shares same column number as the response variable.
   irt_factor   =   0
   IF (KEYWORD_SET(irt)) THEN BEGIN
      IF (IMSL_LONG(irt(0)) EQ (n_class_cvt + n_continuous_cvt)) THEN irt_factor = 0 ELSE irt_factor = 1
   END
   IF (m_cvt NE (KEYWORD_SET(IFREQ)+KEYWORD_SET(IFIX)+KEYWORD_SET(ICEN)+$
                 irt_factor+KEYWORD_SET(ILT) + 1)) THEN $
     MESSAGE,  "The second dimension of X is not the correct size."
   IF ((KEYWORD_SET(var_effects)+KEYWORD_SET(indices_effects)) EQ 1) THEN $
     MESSAGE,  "Keywords VAR_EFFECTS and INDICES_EFFECTS must be used together."
   IF (KEYWORD_SET(var_effects)) THEN BEGIN
      size_ve = IMSL_SIZE(var_effects)
      IF (size_ve(0) NE 1) THEN BEGIN
         message, "VAR_EFFECTS must be a 1-D array."
      END
      n_var_effects_cvt = size_ve(1)
      size_ie = IMSL_SIZE(indices_effects)
      IF (size_ie(0) NE 1) THEN BEGIN
         message, "INDICES_EFFECTS must be a 1-D array."
      END
      IF (size_ie(1) NE IMSL_LONG(TOTAL(var_effects))) THEN $
        MESSAGE,  "The number of elements in INDICES_EFFECTS must equal the sum " + $
        "of the elements OF VAR_EFFECTS"
   END   

   IF (KEYWORD_SET(init_est)) THEN BEGIN
      size_ie = IMSL_SIZE(init_est)
      IF (size_ie(0) NE 1) THEN BEGIN
         message, "INIT_EST must be a 1-D array."
      END
      n_coef_input_cvt  =  size_ie(1)
   END
   est_set  =  ARG_PRESENT(est_prob)+KEYWORD_SET(est_nobs)+KEYWORD_SET(est_time)+KEYWORD_SET(est_npt)+KEYWORD_SET(est_delta)
   IF ((est_set NE 0) AND (est_set NE 5)) THEN MESSAGE,  $
     "Keywords EST_TIME, EST_NPT, EST_DELTA and EST_PROB must be used together."
   IF ((est_set EQ 0) AND (ARG_PRESENT(est_xbeta))) THEN $
     MESSAGE, "Keyword EST_XBETA cannot be used unless EST_PROB is also set."
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
   ; Input LONG keyword(s) and arguments
   model_cvt = IMSL_LONG(model(0))
   IF (KEYWORD_SET(ifreq) EQ TRUE) THEN ifreq_cvt  =  IMSL_LONG(ifreq(0))
   IF (KEYWORD_SET(ifix) EQ TRUE) THEN ifix_cvt  =  IMSL_LONG(ifix(0))
   IF (KEYWORD_SET(icen) EQ TRUE) THEN icen_cvt  =  IMSL_LONG(icen(0))
   IF (KEYWORD_SET(irt) EQ TRUE) THEN irt_cvt  =  IMSL_LONG(irt(0))
   IF (KEYWORD_SET(ilt) EQ TRUE) THEN ilt_cvt  =  IMSL_LONG(ilt(0))
   IF (ARG_PRESENT(itmax)) THEN BEGIN
      IF (SIZE(itmax, /type) EQ 0) THEN itmax_cvt = IMSL_LONG(30) ELSE  itmax_cvt  =  IMSL_LONG(itmax(0))
   END ELSE  itmax_cvt = IMSL_LONG(30)
   IF (KEYWORD_SET(lp_max) EQ TRUE) THEN lp_max_cvt  =  IMSL_LONG(lp_max(0))
   IF (KEYWORD_SET(max_class) EQ TRUE) THEN max_class_cvt  =  IMSL_LONG(max_class(0))
   IF (KEYWORD_SET(var_effects) EQ TRUE) THEN var_effects_cvt  =  IMSL_LONG(var_effects)
   IF (KEYWORD_SET(indices_effects) EQ TRUE) THEN i_effects_cvt  =  IMSL_LONG(indices_effects)
   IF (KEYWORD_SET(no_intercept) EQ TRUE) THEN no_intercept_cvt  =  IMSL_1
   IF (KEYWORD_SET(est_npt) EQ TRUE) THEN BEGIN
      est_npt_cvt     =     IMSL_LONG(est_npt(0))
      IF (est_npt_cvt LT 1) THEN MESSAGE,   "Keyword EST_NPT must be positive."
      est_nobs_cvt     =     IMSL_LONG(est_nobs(0))
      IF (est_nobs_cvt LT 1) THEN MESSAGE,   "Keyword EST_NOBS must be positive."
   END
   
   
   ; Output LONG keyword(s)
   IF ((ARG_PRESENT(n_class_vals) EQ TRUE) OR (ARG_PRESENT(n_class_vals) EQ TRUE)) $
     THEN ncv_spc   =   IMSL_LONARR(n_class_cvt)
   IF (ARG_PRESENT(obs_status) EQ TRUE) THEN obs_status_spc = IMSL_LONARR(nobs_cvt)
   IF (ARG_PRESENT(nmissing) EQ TRUE) THEN nmissing_spc = IMSL_LONG(0)
   IF (ARG_PRESENT(iterations) EQ TRUE) THEN niters_spc = IMSL_LONG(0)
   RESULT = IMSL_0
   ;
   ; Need to get NCOEF in order to dimension some output arrays.
   tmp_coefs  =  IMSL_REGRESSORS(x(*,  0:n_class_cvt+n_continuous_cvt-1),  $
                           n_class_cvt,  n_continuous_cvt,  $
                           /DUMMY_METHOD,  var_effects  =  var_effects,  $
                           indices_effects  =  indices_effects)
   tmp_size = IMSL_SIZE(tmp_coefs)

   intcep = IMSL_1
   IF (KEYWORD_SET(no_intercept)) THEN intcep = IMSL_0
   ncoef  =  tmp_size(2) + intcep
   IF (model_cvt NE 0) THEN ncoef  =  ncoef+IMSL_1
   IF (tmp_size(0) EQ 1) THEN BEGIN 
      ncoef  = 1
      nmeans  =  1
   END ELSE BEGIN
      nmeans  =  ncoef-intcep
      IF (model_cvt NE 0) THEN nmeans  =  nmeans-1
   END
   ; Floating point arguments and keywords
   IF (type EQ TYP_DOUBLE) THEN BEGIN
      ; Input
      x_cvt = DOUBLE(TRANSPOSE(x))
      IF (KEYWORD_SET(eps) EQ TRUE) THEN eps_cvt = DOUBLE(eps(0))
      IF (KEYWORD_SET(est_time) EQ TRUE) THEN est_time_cvt = DOUBLE(est_time(0))
      IF (KEYWORD_SET(est_delta) EQ TRUE) THEN est_delta_cvt = DOUBLE(est_delta(0))
      IF (KEYWORD_SET(init_est) EQ TRUE) THEN init_est_cvt = DOUBLE(init_est)
      ; Output
      IF (ARG_PRESENT(class_vals) EQ TRUE) THEN class_vals_spc  =  DBLARR(n_class_cvt*nobs_cvt)
      IF (ARG_PRESENT(coef_stat) EQ TRUE) THEN coef_stat_spc  =  DBLARR(4, ncoef)
      IF (ARG_PRESENT(criterion) EQ TRUE) THEN criterion_spc  =  DOUBLE(0.0)
      IF (ARG_PRESENT(covariances) EQ TRUE) THEN covariances_spc  =  DBLARR(ncoef, ncoef)
      IF (ARG_PRESENT(means) EQ TRUE) THEN means_spc  =  DBLARR(nmeans)
      IF (ARG_PRESENT(case_analysis) EQ TRUE) THEN ca_spc  =  DBLARR(5, nobs_cvt)
      IF (ARG_PRESENT(iterations) EQ TRUE) THEN iterations_spc  =  DBLARR(5, itmax_cvt)
      IF (ARG_PRESENT(last_step) EQ TRUE) THEN last_step_spc  =  DBLARR(ncoef)      
      IF (ARG_PRESENT(est_prob) EQ TRUE) THEN est_prob_spc  =  DBLARR(2*est_nobs_cvt+1, est_npt_cvt)      
      IF (ARG_PRESENT(est_prob) EQ TRUE) THEN est_xbeta_spc  =  DBLARR(est_nobs_cvt)      
   END ELSE BEGIN
      ; Input
      x_cvt = FLOAT(TRANSPOSE(x))
      IF (KEYWORD_SET(eps) EQ TRUE) THEN eps_cvt = FLOAT(eps(0))
      IF (KEYWORD_SET(est_time) EQ TRUE) THEN est_time_cvt = FLOAT(est_time(0))
      IF (KEYWORD_SET(est_delta) EQ TRUE) THEN est_delta_cvt = FLOAT(est_delta(0))
      IF (KEYWORD_SET(init_est) EQ TRUE) THEN init_est_cvt = FLOAT(init_est)
      ; Output
      IF (ARG_PRESENT(class_vals) EQ TRUE) THEN class_vals_spc  =  FLTARR(n_class_cvt*nobs_cvt)
      IF (ARG_PRESENT(coef_stat) EQ TRUE) THEN coef_stat_spc  =  FLTARR(4, ncoef)
      IF (ARG_PRESENT(criterion) EQ TRUE) THEN criterion_spc  =  FLOAT(0.0)
      IF (ARG_PRESENT(covariances) EQ TRUE) THEN covariances_spc  =  FLTARR(ncoef, ncoef)
      IF (ARG_PRESENT(means) EQ TRUE) THEN means_spc  =  FLTARR(nmeans)
      IF (ARG_PRESENT(case_analysis) EQ TRUE) THEN ca_spc  =  FLTARR(5, nobs_cvt)
      IF (ARG_PRESENT(iterations) EQ TRUE) THEN iterations_spc  =  FLTARR(5, itmax_cvt)
      IF (ARG_PRESENT(last_step) EQ TRUE) THEN last_step_spc  =  FLTARR(ncoef)
      IF (ARG_PRESENT(est_prob) EQ TRUE) THEN est_prob_spc  =  FLTARR(2*est_nobs_cvt+1, est_npt_cvt)      
      IF (ARG_PRESENT(est_prob) EQ TRUE) THEN est_xbeta_spc  =  FLTARR(est_nobs_cvt)      
   END
   ;
   ; Call the system function.
   ;
   err_status = 0L
   MATHSTAT_284,type,err_status, $
                    x_cvt, $
                    ca_spc, $
                    class_vals_spc, $
                    coef_stat_spc, $
                    covariances_spc, $
                    criterion_spc, $
                    eps_cvt, $
                    i_effects_cvt, $
                    ifix_cvt, $
                    ifreq_cvt, $
                    init_est_cvt, $
                    ipar_cvt, $
                    itmax_cvt, $
                    last_step_spc, $
                    m_cvt, $
                    max_class_cvt, $
                    means_spc, $
                    model_cvt, $
                    n_class_cvt, $
                    n_coef_input_cvt, $
                    n_continuous_cvt, $
                    n_var_effects_cvt, $
                    ncv_spc, $
                    nmissing_spc, $
                    no_intercept_cvt, $
                    nobs_cvt, $
                    obs_status_spc, $
                    var_effects_cvt, $
                    x_col_dim, $
                    icen_cvt, $
                    irt_cvt, $
                    ilt_cvt, $
                    lp_max_cvt, $
                    niters_spc, $
                    iterations_spc, $
                    est_nobs_cvt, $
                    est_time_cvt, $
                    est_npt_cvt, $
                    est_delta_cvt, $
                    est_prob_spc, $
                    est_xbeta_spc, $
                    result
    ;
   ; Now copy over all output keywords results.
   ;
   IF (ARG_PRESENT(n_class_vals) EQ TRUE) THEN n_class_vals = ncv_spc 
   IF (ARG_PRESENT(obs_status) EQ TRUE) THEN obs_status = obs_status_spc 
   IF (ARG_PRESENT(nmissing) EQ TRUE) THEN nmissing = nmissing_spc 
   IF (ARG_PRESENT(class_vals) EQ TRUE) THEN class_vals = class_vals_spc(0:TOTAL(n_class_vals)-1) 
   IF (ARG_PRESENT(coef_stat) EQ TRUE) THEN coef_stat = TRANSPOSE(coef_stat_spc)  
   IF (ARG_PRESENT(criterion) EQ TRUE) THEN criterion = criterion_spc
   IF (ARG_PRESENT(covariances) EQ TRUE) THEN covariances = covariances_spc
   IF (ARG_PRESENT(means) EQ TRUE) THEN means = means_spc  
   IF (ARG_PRESENT(case_analysis) EQ TRUE) THEN case_analysis = TRANSPOSE(ca_spc)  
   IF (ARG_PRESENT(iterations) EQ TRUE) THEN iterations = TRANSPOSE(iterations_spc(*, 0:niters_spc-1))  
   IF (ARG_PRESENT(last_step) EQ TRUE) THEN last_step = last_step_spc  
   IF (ARG_PRESENT(est_prob) EQ TRUE) THEN est_prob = TRANSPOSE(est_prob_spc)
   IF (ARG_PRESENT(est_xbeta) EQ TRUE) THEN est_xbeta = est_xbeta_spc
   ;
   ; Return.
   ;
   RETURN, result
END

                   
                   
                   

  
      

  
