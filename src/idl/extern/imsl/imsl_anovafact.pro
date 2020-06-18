; $Id: //depot/idl/IDL_70/idl_src/libs/imsl/imsl_6.0/lib/imsl_anovafact.pro#1 $
;
; Copyright (c) 1970-2006, VISUAL NUMERICS Inc. All Rights Reserved.
; This software is confidential information which is proprietary to and a
; trade secret of Visual Numerics Inc.  Use, duplication or disclosure is
; subject to the terms of an appropriate license agreement.
;
FUNCTION imsl_L_binom, n, m
  binom_v = 1.0
  k = m < (n-m)
  xnp1 = double(n+1)
  xkp1 = double(k+1)
  xi = 0.0;
  FOR i = 1, k DO BEGIN
     xi = xi + 1.0
     binom_v = binom_v * (xnp1 - xi) / xi
  END
  RETURN,binom_v
END

FUNCTION imsl_anovafact, n_levels, $                  ;INPUT 1-D array: LONG
                 y, $                            ;INPUT 1-D array: floating point
                 order=order, $                  ;INPUT Scalar LONG
                 double=double, $                ;INPUT Scalar ON/OFF flag
                 pool_inter=pool_inter, $        ;INPUT Scalar ON/OFF flag
                 pure_error=pure_error, $        ;INPUT Scalar ON/OFF flag
                 anova_table=anova_table, $      ;OUTPUT 1-D array: floating point
                 means=means, $                  ;OUTPUT 1-D array: floating point
                 test_effects=test_effects       ;OUTPUT 1-D array: floating point

@imsl_init.pro
   ON_ERROR, on_err_action
   ;
   ; Error checking.
   ;
   ; - Make sure the input arguments are simple arrays.
   ; - Check the length of the second argument passed to this function.
   ;   It should be equal to the product of the elements of the first
   ;   argument.
   ; - PURE_ERROR and POOL_INTER are mutually exclusive.
   ; -  Check that a valid value is given in ORDER.
   ;
   nargs = n_params()
   IF (nargs NE 2) THEN $
         message, "Incorrect number of arguments."
   size_n_levels = IMSL_SIZE(n_levels)
   IF (size_n_levels(0) NE 1) THEN BEGIN
      message, "N_LEVELS must be a 1-D array."
   END
   n_subscripts = IMSL_N_ELEMENTS(n_levels)
   size_y = IMSL_SIZE(y)
   IF (size_y(0) NE 1) THEN BEGIN
      message, "Y must be a 1-D array."
   END
   exp_lgth_y = IMSL_1
   FOR i = IMSL_0, N_elements(n_levels)-1 DO exp_lgth_y = exp_lgth_y * (IMSL_LONG(n_levels(i)))
   IF (N_ELEMENTS(y) NE exp_lgth_y) THEN $
     message, "The length of Y is incorrect"
   ;
   IF ((KEYWORD_SET(pure_error) EQ TRUE) AND (KEYWORD_SET(pool_inter) EQ TRUE)) THEN $
     message, "PURE_ERROR and POOL_INTER are mutually exclusive."
   IF (KEYWORD_SET(order) EQ true) THEN BEGIN
      IF (((IMSL_LONG(order))(0) LT 1) OR ((IMSL_LONG(order))(0) GT n_subscripts -1)) THEN $
        message, "The value input for ORDER must be in the interval [1, N_ELEMENTS(N)]-1], inclusive" ;
   END
   ;
   ; Decide on what precision to use.
   ; Use the highest precision of the input argument(s).
   ;
   type = TYP_FLOAT
   IF (size_y(N_ELEMENTS(size_y)-2) GT TYP_FLOAT) THEN type = TYP_DOUBLE
   IF (KEYWORD_SET(double) EQ true) THEN type = TYP_DOUBLE
   ;
   ; Compute some stuff that may be needed later on to dimension some
   ; of the output arrays.
   ;
   IF (ARG_PRESENT(means) EQ TRUE) THEN BEGIN
     length_of_means = IMSL_1       ;
     IF (KEYWORD_SET(pool_inter) EQ TRUE) THEN BEGIN
         FOR i = 0, n_subscripts-1 DO BEGIN
            tmp_int = IMSL_LONG(n_levels(i))
            length_of_means = length_of_means*(tmp_int+1) ;
         END
      END ELSE BEGIN
         FOR i = 0, n_subscripts-2 DO BEGIN
            tmp_int = IMSL_LONG(n_levels(i))
            length_of_means = length_of_means*(tmp_int + 1) ;
         END
      END
   END
   ;
   ; Compute the length of TEST_EFFECTS if needed.
   ;
   IF (ARG_PRESENT(test_effects) EQ TRUE) THEN BEGIN
      nef = IMSL_0                  ;
      IF (KEYWORD_SET(pool_inter) EQ TRUE) THEN $
        ntmp = n_subscripts ELSE ntmp = n_subscripts-1
      IF (KEYWORD_SET(order) EQ TRUE) THEN BEGIN
         IF ((IMSL_LONG(order))(0) LT 0) THEN iend = -1*(IMSL_LONG(order))(0) $
             ELSE iend = (IMSL_LONG(order))(0)
      END ELSE iend = n_subscripts-1
      FOR i = 1, iend DO nef = nef + IMSL_LONG(imsl_l_binom(ntmp, i))
   END
   ;
   ; Setup the parameters for the call to the system function.
   ;
   ; Input LONG argument(s) & keyword(s)
   n_levels_cvt = IMSL_LONG(n_levels)
   IF (KEYWORD_SET(order) EQ TRUE) THEN $
     order_cvt = (IMSL_LONG(order))(0)
   IF (KEYWORD_SET(pool_inter) EQ TRUE) THEN $
     pool_inter_cvt = IMSL_1
   IF (KEYWORD_SET(pure_error) EQ TRUE) THEN $
     pure_error_cvt = IMSL_1
   ;
   ; Floating point arguments and keywords
   IF (type EQ TYP_DOUBLE) THEN BEGIN
      result = double(0.0)
      ;
      ; Input
      y_cvt = double(y)
      IF (KEYWORD_SET(confidence) EQ TRUE) THEN $
        confidence_cvt = (double(confidence))(0)
      ;
      ; Output keywords.
      ;
      IF (ARG_PRESENT(anova_table) EQ TRUE) THEN $
        anova_table_spc = dblarr(15)
      IF (ARG_PRESENT(means) EQ TRUE) THEN $
        means_spc = dblarr(length_of_means)
      IF (ARG_PRESENT(test_effects) EQ TRUE) THEN $
        tst_effects_spc = dblarr(4, nef)
   END ELSE BEGIN
      result = float(0.0)
      ;
      ; Input
      y_cvt = float(y)
      IF (KEYWORD_SET(confidence) EQ TRUE) THEN $
        confidence_cvt = (float(confidence))(0)
      ;
      ; Output keywords.
      ;
      IF (ARG_PRESENT(anova_table) EQ TRUE) THEN $
        anova_table_spc = fltarr(15)
      IF (ARG_PRESENT(means) EQ TRUE) THEN $
        means_spc = fltarr(length_of_means)
      IF (ARG_PRESENT(test_effects) EQ TRUE) THEN $
        tst_effects_spc = fltarr(4, nef)
   END
   ;
   ; Call the system function.
   ;
   err_status = 0L
   MATHSTAT_102, type,err_status, n_levels_cvt, y_cvt, n_subscripts, $
                          order_cvt, $
                          pool_inter_cvt, $
                          pure_error_cvt, $
                          anova_table_spc, $
                          means_spc, $
                          tst_effects_spc, $
                          result
   ;
   ; Now copy over all output keywords results.
   ;
   IF (ARG_PRESENT(anova_table) EQ TRUE) THEN $
     anova_table = anova_table_spc
   IF (ARG_PRESENT(means) EQ TRUE) THEN $
     means = means_spc
   IF (ARG_PRESENT(test_effects) EQ TRUE) THEN $
     test_effects = transpose(tst_effects_spc)
   ;
   ; Return.
   ;
   RETURN, result
END
