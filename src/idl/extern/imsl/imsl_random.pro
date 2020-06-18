; $Id: //depot/idl/IDL_70/idl_src/libs/imsl/imsl_6.0/lib/imsl_random.pro#1 $
;
; Copyright (c) 1970-2006, VISUAL NUMERICS Inc. All Rights Reserved.
; This software is confidential information which is proprietary to and a
; trade secret of Visual Numerics Inc.  Use, duplication or disclosure is
; subject to the terms of an appropriate license agreement.
;
FUNCTION imsl_random, n, $                           ;INPUT Scalar LONG
                 a=a, $                         ;INPUT Scalar floating point
                 beta=beta, $                   ;INPUT Scalar ON/OFF flag
                 binomial=binomial, $           ;INPUT Scalar ON/OFF flag
                 cauchy=cauchy, $               ;INPUT Scalar ON/OFF flag
                 chi_squared=chi_squared, $     ;INPUT Scalar ON/OFF flag
                 covariances=covariances, $     ;INPUT 2-D array: floating point
                 double=double, $               ;INPUT Scalar ON/OFF flag
                 discrete_unif=discrete_unif, $ ;INPUT Scalar ON/OFF flag
                 exponential=exponential, $     ;INPUT Scalar ON/OFF flag
                 gamma=gamma, $                 ;INPUT Scalar ON/OFF flag
                 geometric=geometric, $         ;INPUT Scalar ON/OFF flag
                 hypergeometric=hypergeometric, $ ;INPUT Scalar ON/OFF flag
                 logarithmic=logarithmic, $     ;INPUT Scalar ON/OFF flag
                 lognormal=lognormal, $         ;INPUT Scalar ON/OFF flag
                 mix_exponential=mix_exponential, $ ;INPUT Scalar ON/OFF flag
                 mvar_normal=mvar_normal, $     ;INPUT Scalar ON/OFF flag
                 neg_binomial=neg_binomial, $   ;INPUT Scalar ON/OFF flag
                 normal=normal, $               ;INPUT Scalar ON/OFF flag
                 parameters=parameters, $       ;INPUT parameters array (numeric)
                 pin=pin, $                     ;INPUT Scalar floating point
                 poisson=poisson, $             ;INPUT Scalar ON/OFF flag
                 qin=qin, $                     ;INPUT Scalar floating point
                 student_t=student_t, $         ;INPUT Scalar ON/OFF flag
                 theta=theta, $                 ;INPUT Scalar floating point
                 triangular=triangular, $       ;INPUT Scalar ON/OFF flag
                 uniform=uniform, $             ;INPUT Scalar ON/OFF flag
                 von_mises=von_mises, $         ;INPUT Scalar ON/OFF flag
                 weibull=weibull, $             ;INPUT Scalar ON/OFF flag
                 permutation=permutation, $     ;INPUT Scalar ON/OFF flag
                 sphere=sphere, $               ;INPUT Scalar ON/OFF flag
                 stable=stable, $               ;INPUT Scalar ON/OFF flag
                 multinomial=multinomial, $     ;INPUT Scalar ON/OFF flag
                 probabilities=probabilities, $ ;INPUT Scalar ON/OFF flag
                 sample_indices=sample_indices  ;INPUT Scalar ON/OFF flag


@imsl_init.pro
   ON_ERROR, on_err_action
   ;
   ; Error checking:
   ; n must be at least 1.
   ; Figure out which type of random numbers are to be generated.
   ;  dist_type = 1, ..., 21 where
   ;  1 --> Uniform
   ;  2 --> Normal
   ;  3 --> Exponential
   ;  4 --> Poisson
   ;  5 --> Gamma
   ;  6 --> Beta
   ;  7 --> Normal Multivariate
   ;  8 --> Triangular
   ;  9 --> Cauchy
   ;  10 --> chi_squared
   ;  11 --> von Mises
   ;  12 --> Lognormal
   ;  13 --> Weibull
   ;  14 --> Student t
   ;  15 --> Mixed Exponential
   ;  16 --> Discrete uniform
   ;  17 --> Geometric
   ;  18 --> Logarithmic
   ;  19 --> negative binomial
   ;  20 --> hypergeometric
   ;  21 --> binomial
   ;  22 --> permutation
   ;  23 --> sphere
   ;  24 --> sample indices
   ;  25 --> stable
   ;  26 --> multinomial
   ;
   ; Check the keyword associated with dist_type.
   ;
   ;
   nargs = n_params()
   IF (nargs NE 1) THEN $
     message, 'Incorrect number of arguments.'
   n_cvt = IMSL_LONG(n(0))
   IF (n_cvt LT 1) THEN message, 'The number of points must be greater than zero.'

   dist_type_vect = [ KEYWORD_SET(uniform), $
                      KEYWORD_SET(normal), $
                      KEYWORD_SET(exponential), $
                      KEYWORD_SET(poisson), $
                      KEYWORD_SET(gamma), $
                      KEYWORD_SET(beta), $
                      KEYWORD_SET(mvar_normal), $
                      KEYWORD_SET(triangular), $
                      KEYWORD_SET(cauchy), $
                      KEYWORD_SET(chi_squared), $
                      KEYWORD_SET(von_mises), $
                      KEYWORD_SET(lognormal), $
                      KEYWORD_SET(weibull), $
                      KEYWORD_SET(student_t), $
                      KEYWORD_SET(mix_exponential), $
                      KEYWORD_SET(discrete_unif), $
                      KEYWORD_SET(geometric), $
                      KEYWORD_SET(logarithmic), $
                      KEYWORD_SET(neg_binomial), $
                      KEYWORD_SET(hypergeometric), $
                      KEYWORD_SET(binomial), $
                      KEYWORD_SET(permutation), $
                      KEYWORD_SET(sphere), $
                      KEYWORD_SET(sample_indices), $
                      KEYWORD_SET(stable), $
                      KEYWORD_SET(multinomial) ]

   ;
   ;Select the type of random numbers to produce, based upon the keywords.
   ;Default is uniform distribution
   ;Only one 'type' keyword can be supplied at a time
   CASE (TOTAL(dist_type_vect)) OF
      0:    BEGIN & uniform = 1 & dist_type = IMSL_1 & END
      1:    dist_type = IMSL_LONG((0 > (where(dist_type_vect EQ 1))(0)) + 1)
      ELSE: MESSAGE, 'Conflicting keywords have been supplied.'
   ENDCASE

   ;
   ;Select the precision to use.

   IF (KEYWORD_SET(double)) THEN $
      calc_type = TYP_DOUBLE $
   ELSE $
      calc_type = TYP_FLOAT

   result_type = calc_type

   ;Test for presence of required keywords.
   CASE (dist_type) OF
      2: BEGIN                  ;===NORMAL===
         IF (KEYWORD_SET(parameters)) THEN nar_cvt = IMSL_1
      END
      4: BEGIN                  ;===Poisson===
         IF (N_ELEMENTS(theta) NE 0) THEN $
           parameters = theta(0)
         IF (N_ELEMENTS(parameters) EQ 0) THEN $
           MESSAGE, "POISSON keyword requires PARAMETERS keyword"

         theta_cvt  = FLOAT(parameters(0))
         calc_type   = TYP_MEMINT
         result_type = TYP_MEMINT
      END
      5: BEGIN                  ;===Gamma===
         IF (N_ELEMENTS(a) NE 0) THEN $
           parameters = a
         If (N_ELEMENTS(parameters) EQ 0) THEN $
           MESSAGE, "GAMMA keyword requires PARAMETERS keyword"
         CASE calc_type OF
           TYP_FLOAT:  a_cvt = FLOAT(parameters(0))
           TYP_DOUBLE: a_cvt = DOUBLE(parameters(0))
         ENDCASE
      END
      6: BEGIN                  ;===Beta===
         IF (N_ELEMENTS(pin) NE 0) THEN BEGIN
           CASE (N_ELEMENTS(parameters)) OF
             0:    parameters = pin(0)
             ELSE: parameters(0) = pin(0)
           ENDCASE
         ENDIF
         IF (N_ELEMENTS(qin) NE 0) THEN BEGIN
           CASE (N_ELEMENTS(parameters)) OF
             0:
             1:    parameters    = [parameters(0), qin(0)]
             ELSE: parameters(1) = qin(0)
           ENDCASE
         ENDIF
         IF (N_ELEMENTS(parameters) LT 2) THEN $
           MESSAGE, "BETA keyword requires PARAMETERS keyword"
         CASE calc_type OF
           TYP_FLOAT:  BEGIN
               pin_cvt = FLOAT(parameters(0))
               qin_cvt = FLOAT(parameters(1))
             END
           TYP_DOUBLE: BEGIN
               pin_cvt = DOUBLE(parameters(0))
               qin_cvt = DOUBLE(parameters(1))
             END
         ENDCASE
      END
      7: BEGIN                  ;===Normal Multivariate===
         IF (N_ELEMENTS(covariances) EQ 0) THEN $
           MESSAGE, "MVAR_NORMAL keyword requires the COVARIANCES keyword"
         ;Check that COVARIANCES is a 2-D square array.
         size_cov = IMSL_SIZE(covariances)
         IF (size_cov(0) NE 2) THEN $
           MESSAGE, "Covariance array must be 2-D"
         IF (size_cov(1) NE size_cov(2)) THEN $
           MESSAGE, "Covariance array must be square"
         cov_dim = IMSL_LONG(size_cov(1))
         CASE calc_type OF
           TYP_FLOAT:  cov_cvt = FLOAT(covariances)
           TYP_DOUBLE: cov_cvt = DOUBLE(covariances)
         ENDCASE
      END
      10: BEGIN                  ;===Chi Squared===
         IF (N_ELEMENTS(parameters) EQ 0) THEN $
           MESSAGE, "CHI_SQUARED keyword requires the PARAMETERS keyword"
         CASE calc_type OF
           TYP_FLOAT:  df_cvt = FLOAT(parameters(0))
           TYP_DOUBLE: df_cvt = DOUBLE(parameters(0))
         ENDCASE
      END
      11: BEGIN                  ;===von Mises===
         IF (N_ELEMENTS(parameters) EQ 0) THEN $
           MESSAGE, "VON_MISES keyword requires the PARAMETERS keyword"
         c_vm_cvt = parameters(0)
         CASE calc_type OF
           TYP_FLOAT:  c_vm_cvt = FLOAT(parameters(0))
           TYP_DOUBLE: c_vm_cvt = DOUBLE(parameters(0))
         ENDCASE
      END
      12: BEGIN                 ;===Lognormal===
         IF (N_ELEMENTS(parameters) LT 2) THEN $
           MESSAGE, "LOGNORMAL keyword requires 2 PARAMETERS"
         CASE calc_type OF
           TYP_FLOAT:  BEGIN
               mean_log_cvt = FLOAT(parameters(0))
               std_log_cvt  = FLOAT(parameters(1))
             END
           TYP_DOUBLE: BEGIN
               mean_log_cvt = DOUBLE(parameters(0))
               std_log_cvt  = DOUBLE(parameters(1))
             END
         ENDCASE
      END
      13: BEGIN                  ;===Weibull===
         CASE (N_ELEMENTS(parameters)) OF
           0: MESSAGE, "WEIBULL keyword requires the PARAMETERS keyword"
           1: parameters = [DOUBLE(parameters(0)), 1.0D]
           ELSE:
         ENDCASE
         CASE calc_type OF
           TYP_FLOAT:  BEGIN
               a_cvt      = FLOAT(parameters(0))
               b_weib_cvt = FLOAT(parameters(1))
             END
           TYP_DOUBLE: BEGIN
               a_cvt      = DOUBLE(parameters(0))
               b_weib_cvt = DOUBLE(parameters(1))
             END
         ENDCASE
      END
      14: BEGIN                  ;===Student t===
         IF (N_ELEMENTS(parameters) EQ 0) THEN $
           MESSAGE, "STUDENT_T keyword requires the PARAMETERS keyword"
         CASE calc_type OF
           TYP_FLOAT:  df_cvt = FLOAT(parameters(0))
           TYP_DOUBLE: df_cvt = DOUBLE(parameters(0))
         ENDCASE
      END
      15: BEGIN                 ;===Mixed Exponential===
         IF (N_ELEMENTS(parameters) LT 3) THEN $
           MESSAGE, "MIX_EXPONENTIAL keyword requires 3 PARAMETERS"
         CASE calc_type OF
           TYP_FLOAT:  BEGIN
               exp_t1_cvt = FLOAT(parameters(0))
               exp_t2_cvt = FLOAT(parameters(1))
               mix_p_cvt  = FLOAT(parameters(2))
             END
           TYP_DOUBLE: BEGIN
               exp_t1_cvt = DOUBLE(parameters(0))
               exp_t2_cvt = DOUBLE(parameters(1))
               mix_p_cvt  = DOUBLE(parameters(2))
             END
         ENDCASE
      END
      16: BEGIN                  ;===Discrete uniform===
         IF (N_ELEMENTS(parameters) EQ 0) THEN $
           MESSAGE, "DISCRETE_UNIF keyword requires the PARAMETERS keyword"
         disc_k_cvt  = IMSL_LONG(parameters(0))
         result_type = TYP_MEMINT
      END
      17: BEGIN                  ;===Geometric===
         IF (N_ELEMENTS(parameters) EQ 0) THEN $
           MESSAGE, "GEOMETRIC keyword requires the PARAMETERS keyword"
         CASE calc_type OF
           TYP_FLOAT:  p_geom_cvt = FLOAT(parameters(0))
           TYP_DOUBLE: p_geom_cvt = DOUBLE(parameters(0))
         ENDCASE
         result_type = TYP_MEMINT
      END
      18: BEGIN                  ;===Logarithmic===
         IF (N_ELEMENTS(parameters) EQ 0) THEN $
           MESSAGE, "LOGARITHMIC keyword requires the PARAMETERS keyword"
         CASE calc_type OF
           TYP_FLOAT:  a_cvt = FLOAT(parameters(0))
           TYP_DOUBLE: a_cvt = DOUBLE(parameters(0))
         ENDCASE
         result_type = TYP_MEMINT
      END
      19: BEGIN                  ;===Negative Binomial===
         IF (N_ELEMENTS(parameters) LT 2) THEN $
           MESSAGE, "NEG_BINOMIAL keyword requires 2 PARAMETERS"
         nb_rk_cvt = parameters(0)
         nb_p_cvt  = parameters(1)
         CASE calc_type OF
           TYP_FLOAT:  BEGIN
               nb_rk_cvt = FLOAT(parameters(0))
               nb_p_cvt  = FLOAT(parameters(1))
             END
           TYP_DOUBLE: BEGIN
               nb_rk_cvt = DOUBLE(parameters(0))
               nb_p_cvt = DOUBLE(parameters(1))
             END
         ENDCASE
         result_type = TYP_MEMINT
      END
      20: BEGIN                  ;===Hypergeometric===
         IF (N_ELEMENTS(parameters) LT 3) THEN $
           message, "HYPERGEOMETRIC keyword requires 3 PARAMETERS"
         hypergeo_m_cvt = IMSL_LONG(parameters(0))
         hypergeo_n_cvt = IMSL_LONG(parameters(1))
         hypergeo_l_cvt = IMSL_LONG(parameters(2))
         result_type = TYP_MEMINT
      END
      21: BEGIN                  ;===Binomial===
         IF (N_ELEMENTS(parameters) LT 2) THEN $
           MESSAGE, "BINOMIAL keyword requires 2 PARAMETERS"
         CASE calc_type OF
           TYP_FLOAT:  binom_p_cvt = FLOAT(parameters(0))
           TYP_DOUBLE: binom_p_cvt = DOUBLE(parameters(0))
         ENDCASE
         binom_n_cvt = IMSL_LONG(parameters(1))
         result_type = TYP_MEMINT
      END
      22: BEGIN                  ;===Permutation===
         result_type = TYP_MEMINT
      END
      23: BEGIN                  ;===Sphere===
         IF (N_ELEMENTS(parameters) LT 1) THEN $
           MESSAGE, "SPHERE keyword requires 1 PARAMETER"
         kdim_cvt = IMSL_LONG(parameters(0))
      END
      24: BEGIN                  ;===Sample Indices===
         result_type = TYP_MEMINT
         IF (N_ELEMENTS(parameters) LT 1) THEN $
           MESSAGE, "SAMPLE_INDICES keyword requires 1 PARAMETER"
         npopulation_cvt = IMSL_LONG(parameters(0))
      END
      25: BEGIN                  ;===Stable===
         IF (N_ELEMENTS(parameters) LT 2) THEN $
           MESSAGE, "STABLE keyword requires 2 PARAMETERS"
         CASE calc_type OF
           TYP_FLOAT:  BEGIN
               alpha_cvt = FLOAT(parameters(0))
               bprime_cvt  = FLOAT(parameters(1))
             END
           TYP_DOUBLE: BEGIN
               alpha_cvt = DOUBLE(parameters(0))
               bprime_cvt = DOUBLE(parameters(1))
             END
         ENDCASE
	 END
      26: BEGIN                  ;===Multinomial===
         result_type = TYP_MEMINT
         IF (N_ELEMENTS(parameters) LT 1) THEN $
           MESSAGE, "MULTINOMIAL keyword requires 1 PARAMETER"
         ntrials_cvt = IMSL_LONG(parameters(0))
	 if (NOT KEYWORD_SET(probabilities)) THEN $
	   MESSAGE, 'PROBABILITIES must be set if MULTINOMIAL is set."
	 k_cvt = IMSL_N_ELEMENTS(probabilities)
	 if (k_cvt LT 2) THEN $
	   MESSAGE, 'PROBABILITIES must have atleast two elements."
	 probs_cvt = FLOAT(probabilities)
         result_type = TYP_MEMINT
      END
      ELSE:
   END
   ;
   ; Get space for the result in the cases that we return a 1-D vector
   ; of random numbers.
   CASE result_type OF
      TYP_MEMINT:    result = IMSL_LONARR(n_cvt)
      TYP_FLOAT:   result = FLTARR(n_cvt)
      TYP_DOUBLE:  result = DBLARR(n_cvt)
   ENDCASE

   ;
   ; Special case:  Normal Multivariate returns a 2-D array
   IF (dist_type EQ 7) THEN BEGIN
      result = MAKE_ARRAY( size_cov(1), $
                           n_cvt,       $
                           Type=IMSL_SIZE(result(0), /Type))
   ENDIF
   ;
   ; Special case 2:  Sphere returns a 2-D array
   IF (dist_type EQ 23) THEN BEGIN
      result = MAKE_ARRAY(kdim_cvt, $
                           n_cvt,       $
                           Type=IMSL_SIZE(result(0), /Type))
   ENDIF
   ;
   ; Special case 3:  Multinomial returns a 2-D array
   IF (dist_type EQ 26) THEN BEGIN
      result = MAKE_ARRAY(k_cvt, $
                           n_cvt,       $
                           Type=IMSL_SIZE(result(0), /Type))
   ENDIF
   ;
   ; Call the system function.
   ;
   err_status = 0L
   MATHSTAT_181, calc_type, err_status, n_cvt, $
                 a_cvt, $
                 cov_cvt, $
                 cov_dim, $
                 pin_cvt, $
                 qin_cvt, $
                 theta_cvt, $
                 result, $
                 df_cvt, $
                 c_vm_cvt, $
                 mean_log_cvt, $
                 std_log_cvt, $
                 b_weib_cvt, $
                 exp_t1_cvt, $
                 exp_t2_cvt, $
                 mix_p_cvt, $
                 disc_k_cvt, $
                 p_geom_cvt, $
                 nb_rk_cvt, $
                 nb_p_cvt, $
                 hypergeo_m_cvt, $
                 hypergeo_n_cvt, $
                 hypergeo_l_cvt, $
                 binom_p_cvt, $
                 binom_n_cvt, $
                 nar_cvt, $
		 kdim_cvt, $
                 npopulation_cvt, $
		 alpha_cvt, $
		 bprime_cvt, $
                 probs_cvt, $
                 ntrials_cvt, $
                 k_cvt, $
                 dist_type
   ; return
   IF ((dist_type EQ 7) OR (dist_type EQ 23) OR (dist_type EQ 26)) THEN result = transpose(result)
   RETURN, result
END
