; $Id: //depot/idl/IDL_70/idl_src/libs/imsl/imsl_6.0/lib/imsl_laplace_inv.pro#1 $
;
; Copyright (c) 1970-2006, VISUAL NUMERICS Inc. All Rights Reserved.
; This software is confidential information which is proprietary to and a
; trade secret of Visual Numerics Inc.  Use, duplication or disclosure is
; subject to the terms of an appropriate license agreement.
;
FUNCTION imsl_laplace_inv, f, $                 ;INPUT Scalar STRING
                      sigma0,  $           ;INPUT Scalar floating point
                      t,  $                ;INPUT 1-D array floating point
                      pseudo_acc=pseudo_acc, $ ;INPUT Scalar floating point
                      sigma=sigma, $       ;INPUT Scalar floating point
                      bvalue=bvalue, $     ;INPUT Scalar floating point
                      mtop=mtop, $         ;INPUT Scalar LONG
                      err_est=err_est, $   ;OUTPUT Scalar floating point
                      trunc_err=trunc_err, $ ;OUTPUT Scalar floating point
                      cond_err=cond_err, $ ;OUTPUT Scalar floating point
                      disc_err=disc_err, $ ;OUTPUT Scalar floating point
                      k=k, $               ;OUTPUT Scalar floating point
                      r=r, $               ;OUTPUT Scalar floating point
                      big_coef_log=big_coef_log, $ ;OUTPUT Scalar floating point
                      small_coef_log=small_coef_log, $;OUTPUT Scalar floating point
                      indicators=indicators, $ ;OUTPUT 1D array of LONG.
                      double = DOUBLE      ;INPUT Scalar ON/OFF flag

@imsl_init.pro  
   ON_ERROR, on_err_action
   ;
   ; Error checking.
   ; - F must be a scalar string.
   ; - T must be a 1D array, or a scalar.
   ;                       
   nargs = n_params()
   IF (nargs NE 3) THEN message, 'Incorrect number of arguments.'
   size_f = IMSL_SIZE(f)
   IF ((N_ELEMENTS(f) NE 1) OR (size_f(N_ELEMENTS(size_f)-2) NE 7)) THEN $
     message, 'F must be a scalar string.'
   size_t = IMSL_SIZE(t)
   IF (size_t(0) gt 1) THEN message, 'T must be a 1-D array, or a scalar.'
   n = IMSL_N_ELEMENTS(t)

   ; 
   ; Decide on what precision to use.
   type  =  TYP_FLOAT
   size_sigma0= IMSL_SIZE(sigma0)
   IF (size_sigma0(N_ELEMENTS(size_sigma0)-2) EQ  TYP_DOUBLE) THEN type = TYP_DOUBLE
   IF (size_t(N_ELEMENTS(size_t)-2) EQ  TYP_DOUBLE) THEN type = TYP_DOUBLE
   IF (KEYWORD_SET(double) EQ true) THEN type = TYP_DOUBLE
   ;
   ; Setup the parameters for the call to the system function.
   ; Input LONG keyword(s)
   ;
   IF (KEYWORD_SET(mtop)) THEN mtop_cvt = IMSL_LONG(mtop(0)) ELSE mtop_cvt = IMSL_LONG(1024)
   indicators_spc = IMSL_LONARR(n)
   ;
   ; Floating point arguments and keywords   
   IF (type EQ TYP_DOUBLE) THEN BEGIN
      tmp = imsl_machine(/double)
      result = DBLARR(n)
      sigma0_cvt = double(sigma0(0))
      t_cvt = double(t)
      IF (KEYWORD_SET(pseudo_acc)) THEN pseudo_acc_cvt = DOUBLE(pseudo_acc(0)) $
        ELSE pseudo_acc_cvt = SQRT(tmp.(3))
      IF (KEYWORD_SET(sigma)) THEN sigma_cvt =  DOUBLE(sigma(0)) $
        ELSE sigma_cvt = sigma0_cvt+.7
      IF (KEYWORD_SET(bvalue)) THEN bvalu_cvt = DOUBLE(bvalue(0)) $
        ELSE bvalue_cvt = DOUBLE(2.5*(sigma_cvt-sigma0_cvt))
      err_est_spc = double(0.0)
      trunc_err_spc = double(0.0)
      cond_err_spc = double(0.0)
      disc_err_spc = double(0.0)
      k_spc = double(0.0)
      r_spc = double(0.0)
      b_coef_log_spc = double(0.0)
      s_coef_log_spc = double(0.0)
   END ELSE BEGIN
      tmp = imsl_machine(/float)
      result = fltARR(n)
      sigma0_cvt = float(sigma0(0))
      t_cvt = float(t)
      IF (KEYWORD_SET(pseudo_acc)) THEN pseudo_acc_cvt = FLOAT(pseudo_acc(0)) $
        ELSE pseudo_acc_cvt = SQRT(tmp.(3))
      IF (KEYWORD_SET(sigma)) THEN sigma_cvt =  FLOAT(sigma(0)) $
        ELSE sigma_cvt = sigma0_cvt+.7
      IF (KEYWORD_SET(bvalue)) THEN bvalu_cvt = FLOAT(bvalue(0)) $
        ELSE bvalue_cvt = FLOAT(2.5*(sigma_cvt-sigma0_cvt))
      err_est_spc = float(0.0)
      trunc_err_spc = float(0.0)
      cond_err_spc = float(0.0)
      disc_err_spc = float(0.0)
      k_spc = float(0.0)
      r_spc = float(0.0)
      b_coef_log_spc = float(0.0)
      s_coef_log_spc = float(0.0)
   END
   ;
   ; Call the system function.
   ;
   err_status = 0L
   MATHSTAT_223,   type,   err_status,   $
                              f, $
                              sigma0_cvt, $
                              n, $
                              t_cvt, $
                              pseudo_acc_cvt, $
                              sigma_cvt, $
                              bvalue_cvt, $
                              mtop_cvt, $
                              indicators_cvt, $
                              err_est_spc, $
                              trunc_err_spc, $
                              cond_err_spc, $
                              disc_err_spc, $
                              k_spc, $
                              r_spc, $
                              b_coef_log_spc, $
                              s_coef_log_spc, $
                              result
   
      indicators = indicators_spc
      err_est = err_est_spc
      trun_err = trunc_err_spc
      cond_err = cond_err_spc
      disc_err = disc_err_spc
      k = k_spc
      r = r_spc
      big_coef_log = b_coef_log_spc
      small_coef_log = s_coef_log_spc
 RETURN, result
END

