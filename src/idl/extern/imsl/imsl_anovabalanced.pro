; $Id: //depot/idl/IDL_70/idl_src/libs/imsl/imsl_6.0/lib/imsl_anovabalanced.pro#1 $
;
; Copyright (c) 1970-2006, VISUAL NUMERICS Inc. All Rights Reserved.
; This software is confidential information which is proprietary to and a
; trade secret of Visual Numerics Inc.  Use, duplication or disclosure is
; subject to the terms of an appropriate license agreement.
;

FUNCTION imsl_anovabalanced, n_levels  , $    ;INPUT 1-D array: LONG
                     y, $                     ;INPUT 1-D array: floating point
                     n_random, $              ;INPUT Scalar LONG
                     idx_rand_fct, $          ;INPUT 1-D array: LONG
                     n_fct_per_eff, $         ;INPUT 1-D array: LONG
                     idx_fct_per_eff, $       ;INPUT 1-D array: LONG
                     double=double, $         ;INPUT Scalar ON/OFF flag
                     anova_table=anova_table, $ ;OUTPUT 1-D array: floating point
                     confidence=confidence, $ ;INPUT 1-D Scalar floating point
                     model=model, $           ;INPUT 1-D Scalar LONG
                     ems=ems, $               ;OUTPUT 1-D array: floating point
                     y_means=y_means, $       ;OUTPUT 1-D array: floating point
                     var_comp=var_comp        ;OUTPUT 2-D array: floating point

@imsl_init.pro
   ON_ERROR, on_err_action

   ;
   ; Error checking:
   ; The following checks are performed.
   ;   N_LEVELS must be a 1D array, N_FACTORS is set to the length.
   ;   Y must be a 1D array.
   ;   The length of Y must be the product of the elements of N_LEVELS.
   ;   IDX_RAND_FAC must be a 1D array of length ABS(N_RANDOM)
   ;   N_FCT_PER_EFF must be a 1D array, N_MODEL_EFFECTS is set to the length.
   ;   IDX_FCT_PER_EFF must be a 1D array.
   ;   The length of IDX_FCT_PER_EFF must be the sum of the elements of N_FCT_PER_EFF.
   ;
   nargs = n_params()
   IF (nargs NE 6) THEN message, "Incorrect number of arguments."

   size_nl = IMSL_SIZE(n_levels)
   IF (size_nl(0) NE 1) THEN BEGIN
      message, "N_LEVELS must be a 1-D array."
   END
   n_factors = size_nl(1)
   pnl  =  IMSL_1
   FOR i  =  0,  n_factors-1 DO pnl  =  pnl*n_levels(i)
   
   size_y = IMSL_SIZE(y)
   IF (size_y(0) NE 1) THEN BEGIN
      message, "Y must be a 1-D array."
   END
   IF (size_y(1) NE pnl) THEN MESSAGE,  "Y is not the correct size."
   
   size_idx_rand = IMSL_SIZE(idx_rand_fct)
   IF (size_idx_rand(0) NE 1) THEN BEGIN
      message, "IDX_RAND_FAC must be a 1-D array."
   END
   IF (size_idx_rand(1) NE ABS(n_random)) THEN MESSAGE,  "IDX_RAND_FCT is not the correct size."

   size_n_fct_per= IMSL_SIZE(n_fct_per_eff)
   IF (size_n_fct_per(0) NE 1) THEN BEGIN
      message, "N_FCT_PER_EFF must be a 1-D array."
   END
   n_model_eff = size_n_fct_per(1)

   size_idx_fct_per= IMSL_SIZE(idx_fct_per_eff)
   IF (size_idx_fct_per(0) NE 1) THEN BEGIN
      message, "IDX_FCT_PER_EFF must be a 1-D array."
   END
   IF (size_idx_fct_per(1) NE TOTAL(n_fct_per_eff)) THEN MESSAGE,  "IDX_FCT_PER_EFF is not the correct size."
   ;
   ; Decide on what precision to use.
   ; Use the highest precision of the input argument(s).
   ;
   type = TYP_FLOAT
   IF (size_y(N_ELEMENTS(size_y)-2) GT TYP_FLOAT) THEN type = TYP_DOUBLE
   IF (KEYWORD_SET(double) EQ true) THEN type = TYP_DOUBLE
   ;
   ; Setup the parameters for the call to the system function.
   ;
   ; Input LONG keyword(s) and arguments
   n_levels_cvt  =  IMSL_LONG(n_levels)
   n_random_cvt  =  IMSL_LONG(n_random(0))
   irf_cvt  =  IMSL_LONG(idx_rand_fct)
   nfpe_cvt  =  IMSL_LONG(n_fct_per_eff)
   ifpe_cvt   =   IMSL_LONG(idx_fct_per_eff)
   IF (KEYWORD_SET(model) EQ TRUE) THEN model_cvt = IMSL_LONG(model(0))
   ; Output LONG keyword(s)
   ;
   ; Compute length of YMEANS.
   lindef = TOTAL(n_fct_per_eff)
   n = n_factors-1
   for i=0, lindef-1 DO BEGIN
      if (idx_fct_per_eff(i) EQ n_factors) THEN BEGIN
	n = n_factors
	i = lindef
      ENDIF
   END
   nymeans  =  IMSL_1
   FOR i  =  0,  n-1 DO nymeans  =  nymeans * (n_levels_cvt(i) + 1)
   
   ; Floating point arguments and keywords
   IF (type EQ TYP_DOUBLE) THEN BEGIN
      ; Input
      y_cvt = double(y)
      IF (KEYWORD_SET(confidence) EQ TRUE) THEN conf_cvt = DOUBLE(confidence) ELSE conf_cvt = DOUBLE(95.)
      ; Output
      IF (ARG_PRESENT(anova_table) EQ TRUE) THEN at_spc = DBLARR(15)
      IF (ARG_PRESENT(var_comp) EQ TRUE) THEN vc_spc = DBLARR(9, n_model_eff+1)
      IF (ARG_PRESENT(ems) EQ TRUE) THEN ems_spc = DBLARR(((n_model_eff+1)*(n_model_eff+2))/2)
      IF (ARG_PRESENT(y_means) EQ TRUE) THEN y_means_spc= DBLARR(nymeans)
      result = DOUBLE(0.0)
   END ELSE BEGIN
      ; Input
      y_cvt = FLOAT(y)
      IF (KEYWORD_SET(confidence) EQ TRUE) THEN conf_cvt = FLOAT(confidence) ELSE conf_cvt = FLOAT(95.)
      ; Output
      IF (ARG_PRESENT(anova_table) EQ TRUE) THEN at_spc = FLTARR(15)
      IF (ARG_PRESENT(var_comp) EQ TRUE) THEN vc_spc = FLTARR(9, n_model_eff+1)
      IF (ARG_PRESENT(ems) EQ TRUE) THEN ems_spc = FLTARR(((n_model_eff+1)*(n_model_eff+2))/2)
      IF (ARG_PRESENT(y_means) EQ TRUE) THEN y_means_spc= FLTARR(nymeans)
      result = FLOAT(0.0)
   END
   ;
   ; Call the system function.
   ;
   err_status = 0L
   MATHSTAT_276,type,err_status, $
                           n_levels_cvt, $
                           n_random_cvt, $
                           irf_cvt, $
                           n_model_eff, $
                           nfpe_cvt, $
                           ifpe_cvt, $
                           n_factors, $
                           model_cvt, $
                           y_cvt, $
                           conf_cvt, $
                           at_spc, $
                           vc_spc, $
                           ems_spc, $
                           y_means_spc, $
                           result
   ;
   ; Now copy over all output keywords results.
   ;
   IF (ARG_PRESENT(anova_table) EQ TRUE) THEN anova_table = at_spc
   IF (ARG_PRESENT(var_comp) EQ TRUE) THEN var_comp = TRANSPOSE(vc_spc)
   IF (ARG_PRESENT(ems) EQ TRUE) THEN ems = ems_spc
   IF (ARG_PRESENT(y_means) EQ TRUE) THEN y_means = y_means_spc
   ;
   ; Return.
   ;
   RETURN, result
END

                   
                   
                   

  
      

  
