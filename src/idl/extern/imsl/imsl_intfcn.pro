; $Id: //depot/idl/IDL_70/idl_src/libs/imsl/imsl_6.0/lib/imsl_intfcn.pro#1 $
;
; Copyright (c) 1970-2006, VISUAL NUMERICS Inc. All Rights Reserved.
; This software is confidential information which is proprietary to and a
; trade secret of Visual Numerics Inc.  Use, duplication or disclosure is
; subject to the terms of an appropriate license agreement.
;
FUNCTION imsl_intfcn, f, $                           ;INPUT Scalar STRING
                   arg1, $
                   arg2, $
                   arg3, $
                   arg4, $
                   algebraic=algebraic, $         ;INPUT Scalar ON/OFF flag
                   alg_left_log=alg_left_log, $   ;INPUT Scalar ON/OFF flag
                   alg_log=alg_log, $             ;INPUT Scalar ON/OFF flag
                   alg_right_log=alg_right_log, $ ;INPUT Scalar ON/OFF flag
                   bound_inf=bound_inf, $         ;INPUT Scalar ON/OFF flag
                   cauchy=cauchy, $               ;INPUT Scalar ON/OFF flag
                   cosine=cosine, $               ;INPUT Scalar ON/OFF flag
                   double=double, $               ;INPUT Scalar ON/OFF flag
                   err_abs=err_abs, $             ;INPUT Scalar floating point
                   err_est=err_est, $             ;OUTPUT Scalar floating point
                   err_rel=err_rel, $             ;INPUT Scalar floating point
                   inf_bound=inf_bound, $         ;INPUT Scalar ON/OFF flag
                   inf_inf=inf_inf, $             ;INPUT Scalar ON/OFF flag
                   max_cycles=max_cycles, $       ;INPUT Scalar ON/OFF flag
                   max_moments=max_moments, $     ;INPUT Scalar ON/OFF flag
                   max_subinter=max_subinter, $   ;INPUT Scalar ON/OFF flag
                   n_cycles=n_cycles, $           ;OUTPUT Scalar LONG
                   n_evals=n_evals, $             ;OUTPUT Scalar LONG
                   n_subinter=n_subinter, $       ;OUTPUT Scalar LONG
                   rule=rule, $                   ;INPUT Scalar LONG
                   sine=sine, $                   ;INPUT Scalar ON/OFF flag
                   sing_pts=sing_pts, $           ;INPUT 1-D array: floating point
                   smooth=smooth, $               ;INPUT Scalar ON/OFF flag
                   two_dimensional=two_dimensional;INPUT Scalar ON/OFF flag

@imsl_init.pro   
   ON_ERROR, on_err_action
   
   max_cycles_cvt = IMSL_LONG(500)
   max_sub_cvt = IMSL_LONG(500)
   max_moments_cvt = IMSL_LONG(21)
   n_subinter_spc = IMSL_0
   n_evals_spc = IMSL_0
   INT_FCN_SING     = 0
   INT_FCN          = 1
   INT_FCN_SING_PTS = 2
   INT_FCN_ALG_LOG  = 3
   INT_FCN_INF      = 4
   INT_FCN_TRIG     = 5
   INT_FCN_FOURIER  = 6
   INT_FCN_CAUCHY   = 7
   INT_FCN_SMOOTH   = 8
   INT_FCN_2D       = 9
   INT_FCN_INCONSISTENT=10  ; implies the arguments and keywords do not imply a valid method.
   ;
   ; Error checking.
   ; This routine can be called with anywhere from one to five
   ; positional arguments.  The first thing to do is decide what
   ; method of quadrature the keywords imply.  Here is the method
   ; used to go through the keywords and decide what method is used.
   ; If:
   ;  o  nargs EQ 1:
   ;     - The keyword INF_INF must be set.
   ;  o  nargs EQ 2:
   ;     - One of the keywords INF_BOUND or BOUND_INF must be set.
   ;  o  nargs EQ 3:
   ;     - (1) If RULE is set, integrate using Gauss-Kronrod rules.
   ;     - (2) If SING_PTS is set, then integrate a function with singular
   ;           points given.
   ;     - (3) If SINE or COSINE is set, then compute the Fourier transform.
   ;     - (4) If SMOOTH is set, then use a nonadaptive method.
   ;     - (5) If none of the above keywords is set, then use the default method.
   ;  o  nargs EQ 4:
   ;     - (1) If CAUCHY is set, the integrate in the Cauch Principle Value sense.
   ;     - (2) If one of SINE or COSINE is set, then integrate a function with a sine
   ;           or cosine factor.
   ;  o  nargs EQ 5:
   ;     - (1) If TWO_DIMENSIONAL is set, then compute a two-dimensional iterated
   ;           integral.
   ;     - (2) If one of ALGEBRAIC, ALG_LEFT_LOG, ALG_LOG, or ALG_RIGHT_LOG is set
   ;           then integrate a function with algebraic singularities.
   ;
   ;  Once a method is decided upon, the relevant keywords ad arguments are checked.
   ;  This means that keywords that do not apply the method that is going to be 
   ;  used are ignored.
   ;                       
   nargs = n_params()
   itgr_method = INT_FCN_INCONSISTENT
   IF (KEYWORD_SET(double)) THEN type = TYP_DOUBLE ELSE type = TYP_FLOAT
   CASE nargs OF
      1: BEGIN
         IF (KEYWORD_SET(inf_inf)) THEN itgr_method = INT_FCN_INF
      END
      2: BEGIN
         size_tmp = IMSL_SIZE(arg1)
         IF (size_tmp(N_ELEMENTS(size_tmp)-2) EQ TYP_DOUBLE) THEN type = TYP_DOUBLE
         IF ((KEYWORD_SET(bound_inf) OR (KEYWORD_SET(inf_bound)))) THEN itgr_method = INT_FCN_INF
      END
      3: BEGIN
         size_tmp = IMSL_SIZE(arg1)
         IF (size_tmp(N_ELEMENTS(size_tmp)-2) EQ TYP_DOUBLE) THEN type = TYP_DOUBLE
         size_tmp = IMSL_SIZE(arg2)
         IF (size_tmp(N_ELEMENTS(size_tmp)-2) EQ TYP_DOUBLE) THEN type = TYP_DOUBLE
         itgr_method = INT_FCN_SING
         IF (KEYWORD_SET(smooth)) THEN itgr_method = INT_FCN_SMOOTH
         IF (KEYWORD_SET(sine) OR KEYWORD_SET(cosine)) THEN itgr_method = INT_FCN_FOURIER
         IF (KEYWORD_SET(sing_pts_set)) THEN itgr_method = INT_FCN_SING_PTS
         IF (KEYWORD_SET(rule)) THEN itgr_method = INT_FCN
      END
      4: BEGIN
         size_tmp = IMSL_SIZE(arg1)
         IF (size_tmp(N_ELEMENTS(size_tmp)-2) EQ TYP_DOUBLE) THEN type = TYP_DOUBLE
         size_tmp = IMSL_SIZE(arg2)
         IF (size_tmp(N_ELEMENTS(size_tmp)-2) EQ TYP_DOUBLE) THEN type = TYP_DOUBLE
         size_tmp = IMSL_SIZE(arg3)
         IF (size_tmp(N_ELEMENTS(size_tmp)-2) EQ TYP_DOUBLE) THEN type = TYP_DOUBLE
         itgr_method = INT_FCN_INCONSISTENT
         IF (KEYWORD_SET(cauchy)) THEN itgr_method = INT_FCN_CAUCHY
         IF (KEYWORD_SET(sine) OR KEYWORD_SET(cosine)) THEN itgr_method = INT_FCN_TRIG
      END
      5: BEGIN
         size_tmp = IMSL_SIZE(arg1)
         IF (size_tmp(N_ELEMENTS(size_tmp)-2) EQ TYP_DOUBLE) THEN type = TYP_DOUBLE
         size_tmp = IMSL_SIZE(arg2)
         IF (size_tmp(N_ELEMENTS(size_tmp)-2) EQ TYP_DOUBLE) THEN type = TYP_DOUBLE
         size_tmp = IMSL_SIZE(arg3)
         IF (size_tmp(N_ELEMENTS(size_tmp)-2) EQ TYP_DOUBLE) THEN type = TYP_DOUBLE
         size_tmp = IMSL_SIZE(arg4)
         IF (size_tmp(N_ELEMENTS(size_tmp)-2) EQ TYP_DOUBLE) THEN type = TYP_DOUBLE
         itgr_method = INT_FCN_INCONSISTENT
         IF (KEYWORD_SET(two_dimensional)) THEN  itgr_method = INT_FCN_2D
         IF ((KEYWORD_SET(algebraic) +  KEYWORD_SET(alg_left_log) + $
              KEYWORD_SET(alg_log) + KEYWORD_SET(alg_right_log)) GT 0) $
           THEN itgr_method = INT_FCN_ALG_LOG
      END
      ELSE: message, 'Incorrect number of arguments.'
   END
   ;
   ; The method is now decided. 
   ; Time to check the relevant arguments and keywords. 
   ; In all cases, p_argv[0] must be a scalar string. 
   ;
   size_f = IMSL_SIZE(f)
   IF ((N_ELEMENTS(f) NE 1) OR (size_f(N_ELEMENTS(size_f)-2) NE 7)) THEN $
     message, 'F must be a scalar string.'

   CASE itgr_method OF
      INT_FCN_SING: BEGIN
         IF (type EQ TYP_DOUBLE) THEN a_cvt = double(arg1(0)) ELSE a_cvt = float(arg1(0))
         IF (type EQ TYP_DOUBLE) THEN b_cvt = double(arg2(0)) ELSE b_cvt = float(arg2(0))
      END
      INT_FCN:BEGIN
         IF (type EQ TYP_DOUBLE) THEN a_cvt = double(arg1(0)) ELSE a_cvt = float(arg1(0))
         IF (type EQ TYP_DOUBLE) THEN b_cvt = double(arg2(0)) ELSE b_cvt = float(arg2(0))
         rule_cvt = IMSL_LONG(rule(0))
      END
      INT_FCN_SING_PTS:BEGIN
         IF (type EQ TYP_DOUBLE) THEN a_cvt = double(arg1(0)) ELSE a_cvt = float(arg1(0))
         IF (type EQ TYP_DOUBLE) THEN b_cvt = double(arg2(0)) ELSE b_cvt = float(arg2(0))
         size_sing_pts = IMSL_SIZE(sing_pts)
         IF (size_sing_pts(0) NE 1) THEN message, 'SING_PTS must be a 1-D array.'
         IF (type EQ TYP_DOUBLE) THEN sing_pts_cvt = double(sing_pts) $
           ELSE sing_pts_cvt = float(sing_pts)
         n_sing_pts = IMSL_N_ELEMENTS(sing_pts)
      END
      INT_FCN_ALG_LOG:BEGIN
         IF (type EQ TYP_DOUBLE) THEN a_cvt = double(arg1(0)) ELSE a_cvt = float(arg1(0))
         IF (type EQ TYP_DOUBLE) THEN b_cvt = double(arg2(0)) ELSE b_cvt = float(arg2(0))
         IF (type EQ TYP_DOUBLE) THEN alpha_cvt = double(arg3(0)) ELSE alpha_cvt = float(arg3(0))
         IF (type EQ TYP_DOUBLE) THEN beta_cvt = double(arg4(0)) ELSE beta_cvt = float(arg4(0))
         IF ((KEYWORD_SET(algebraic) +  KEYWORD_SET(alg_left_log) + $
              KEYWORD_SET(alg_log) + KEYWORD_SET(alg_right_log)) GT 1) $
           THEN message, 'The keywords ALGEBRAIC, ALG_LEFT_LOG, ALG_LOG, AND ' + $
           'ALG_RIGHT_LOG are mutually exclusive.'
      END
      INT_FCN_INF: BEGIN
         IF (type EQ TYP_DOUBLE) THEN BEGIN 
            IF (nargs EQ 2) THEN bound_cvt = double(arg1(0)) ELSE bound_cvt = double(0.0)
         END ELSE BEGIN
            IF (nargs EQ 2) THEN bound_cvt = float(arg1(0)) ELSE bound_cvt = float(0.0)
         END
         IF ((KEYWORD_SET(inf_inf) + KEYWORD_SET(bound_inf) + $
              KEYWORD_SET(inf_bound)) GT 1) THEN message, $
           'The keywords INF_INF, INF_BOUND, and BOUND_INF are mutually exclusive.'
      END
      INT_FCN_TRIG:BEGIN
         IF (type EQ TYP_DOUBLE) THEN a_cvt = double(arg1(0)) ELSE a_cvt = float(arg1(0))
         IF (type EQ TYP_DOUBLE) THEN b_cvt = double(arg2(0)) ELSE b_cvt = float(arg2(0))
         IF (type EQ TYP_DOUBLE) THEN omega_cvt = double(arg3(0)) ELSE omega_cvt = float(arg3(0))
         IF (KEYWORD_SET(max_moments)) THEN max_moments_cvt = IMSL_LONG(max_moments(0))
         IF ((KEYWORD_SET(sine) + KEYWORD_SET(cosine)) GT 1) THEN message, $
           'The keywords SINE and COSINE are mutually exclusive.'
      END
      INT_FCN_FOURIER:BEGIN
         IF (type EQ TYP_DOUBLE) THEN a_cvt = double(arg1(0)) ELSE a_cvt = float(arg1(0))
         IF (type EQ TYP_DOUBLE) THEN omega_cvt = double(arg2(0)) ELSE omega_cvt = float(arg2(0))
         IF (KEYWORD_SET(max_moments)) THEN max_moments_cvt = IMSL_LONG(max_moments(0))
         IF (KEYWORD_SET(max_cycles)) THEN max_cycles_cvt = IMSL_LONG(max_cycles(0))
         n_cycles_spc = IMSL_0
         IF ((KEYWORD_SET(sine) + KEYWORD_SET(cosine)) GT 1) THEN message, $
           'The keywords SINE and COSINE are mutually exclusive.'
         ; Although the documentation doesn't say so, the following
         ; keywords are not allowed for this option. */
         IF (KEYWORD_SET(n_subinter)) THEN message, $
           'N_SUBINTER is not valid for the specified method of integration.'
         IF (KEYWORD_SET(err_rel)) THEN message, $
           'ERR_REL is not valid for the specified method of integration.'
      END
      INT_FCN_CAUCHY:BEGIN
         IF (type EQ TYP_DOUBLE) THEN a_cvt = double(arg1(0)) ELSE a_cvt = float(arg1(0))
         IF (type EQ TYP_DOUBLE) THEN b_cvt = double(arg2(0)) ELSE b_cvt = float(arg2(0))
         IF (type EQ TYP_DOUBLE) THEN c_cvt = double(arg3(0)) ELSE c_cvt = float(arg3(0))
         IF (KEYWORD_SET(n_subinter)) THEN message, $
           'N_SUBINTER is not valid for the specified method of integration.'
         IF (KEYWORD_SET(err_rel)) THEN message, $
           'ERR_REL is not valid for the specified method of integration.'
      END
      INT_FCN_SMOOTH:BEGIN
         IF (type EQ TYP_DOUBLE) THEN a_cvt = double(arg1(0)) ELSE a_cvt = float(arg1(0))
         IF (type EQ TYP_DOUBLE) THEN b_cvt = double(arg2(0)) ELSE b_cvt = float(arg2(0))
         IF (KEYWORD_SET(n_subinter)) THEN message, $
           'N_SUBINTER is not valid for the specified method of integration.'
         IF (KEYWORD_SET(max_subinter)) THEN message, $
           'MAX_SUBINTER is not valid for the specified method of integration.'
         IF (KEYWORD_SET(n_evals)) THEN message, $
           'N_EVALS is not valid for the specified method of integration.'
      END
      INT_FCN_2D:BEGIN
         IF (type EQ TYP_DOUBLE) THEN a_cvt = double(arg1(0)) ELSE a_cvt = float(arg1(0))
         IF (type EQ TYP_DOUBLE) THEN b_cvt = double(arg2(0)) ELSE b_cvt = float(arg2(0))
         size_g = IMSL_SIZE(arg3)
         IF ((N_ELEMENTS(arg3) NE 1) OR (size_g(N_ELEMENTS(size_g)-2) NE 7)) THEN $
           message, 'G must be a scalar string.'
         g = arg3
         size_h = IMSL_SIZE(arg3)
         IF ((N_ELEMENTS(arg4) NE 1) OR (size_h(N_ELEMENTS(size_h)-2) NE 7)) THEN $
           message, 'H must be a scalar string.'
         h = arg4
      END
      INT_FCN_INCONSISTENT:BEGIN
         message, 'The arguments and keywords are not consistent with any ' + $
           'available integration methods.'
      END
   END
   IF (KEYWORD_SET(max_subinter)) THEN max_sub_cvt = IMSL_LONG(max_subinter(0))
   IF (type EQ TYP_DOUBLE) THEN BEGIN
      tmp = imsl_machine(/double)
      err_abs_cvt = sqrt(tmp.(3))
      err_rel_cvt = sqrt(tmp.(3))
      result = double(0.0)
      IF (KEYWORD_SET(err_abs)) THEN err_abs_cvt = double(err_abs(0))
      IF (KEYWORD_SET(err_rel)) THEN err_rel_cvt = double(err_rel(0))
      err_est_spc = double(0.0)
   END ELSE BEGIN
      tmp = imsl_machine(/float)
      err_abs_cvt = sqrt(tmp.(3))
      err_rel_cvt = sqrt(tmp.(3))
      result = float(0.0)
      IF (KEYWORD_SET(err_abs)) THEN err_abs_cvt = float(err_abs(0))
      IF (KEYWORD_SET(err_rel)) THEN err_rel_cvt = float(err_rel(0))
      err_est_spc = float(0.0)
   END
   ;
   ; Call the system function.
   ;
   err_status = 0L
   MATHSTAT_150, type, err_status, $
                              f, $
                              a_cvt, $
                              b_cvt, $
                              c_cvt, $
                              omega_cvt, $
                              bound_cvt, $
                              alpha_cvt, $
                              beta_cvt, $
                              g, $
                              h, $
                              IMSL_LONG(KEYWORD_SET(algebraic)), $
                              IMSL_LONG(KEYWORD_SET(alg_left_log)), $
                              IMSL_LONG(KEYWORD_SET(alg_log)), $
                              IMSL_LONG(KEYWORD_SET(alg_right_log)), $
                              IMSL_LONG(KEYWORD_SET(bound_inf)), $
                              IMSL_LONG(KEYWORD_SET(cauchy)), $
                              IMSL_LONG(KEYWORD_SET(cosine)), $
                              err_abs_cvt, $
                              err_est_spc, $
                              err_rel_cvt, $
                              IMSL_LONG(KEYWORD_SET(inf_bound)), $
                              IMSL_LONG(KEYWORD_SET(inf_inf)), $
                              max_cycles_cvt, $
                              max_moments_cvt, $
                              max_sub_cvt, $
                              n_cycles_spc, $
                              n_evals_spc, $
                              n_subinter_spc, $
                              rule_cvt, $
                              IMSL_LONG(KEYWORD_SET(sine)), $
                              sing_pts_cvt, $
                              n_sing_pts, $
                              IMSL_LONG(KEYWORD_SET(smooth)), $
                              IMSL_LONG(KEYWORD_SET(two_dimensional)), $
                              IMSL_LONG(itgr_method), $
                              result
   IF (ARG_PRESENT(n_evals)) THEN n_evals = n_evals_spc
   IF (ARG_PRESENT(n_subinter)) THEN n_subinter = n_subinter_spc
   IF (ARG_PRESENT(err_est)) THEN err_est = err_est_spc
   IF (ARG_PRESENT(n_cycles)) THEN n_cycles = n_cycles_spc
 RETURN, result
END

