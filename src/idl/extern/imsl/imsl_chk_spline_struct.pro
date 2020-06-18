; $Id: //depot/idl/IDL_70/idl_src/libs/imsl/imsl_6.0/lib/imsl_chk_spline_struct.pro#1 $
;
; Copyright (c) 1970-2006, VISUAL NUMERICS Inc. All Rights Reserved.
; This software is confidential information which is proprietary to and a
; trade secret of Visual Numerics Inc.  Use, duplication or disclosure is
; subject to the terms of an appropriate license agreement.
;
function imsl_Chk_spline_struct, sp, type, struct_is_ppoly
@imsl_init.pro

   ; The following checks are made on the argument that
   ; is supposed to be a spline structure:
   ; - Must be a structure.
   ; - Must have 7 tags
   ; - The 5-th tag determines if the spline is a ppoly or Bspline.
   ; - The type of the 7-th tag determines if the spline
   ;   TYP_FLOAT or TYP_DOUBLE.
   ; - Tagnames must match expected tagnames.
   ; - For each tag:
   ;       o  check the data type.
   ;       o  check the size.
   ;
   err_str = "SPLINE must be a valid spline structure."
   size_sp = IMSL_SIZE(sp)
   IF ((size_sp(0) NE 1) OR $
       (size_sp(N_ELEMENTS(size_sp)-2) NE 8)) THEN $
     message, err_str

   n_tags = n_tags(sp)
   IF (n_tags NE 7) THEN message, err_str
   sp_tagnames = tag_names(sp)
   type = SIZE(sp.(6), /TYPE)
   CASE sp_tagnames(4) OF
      "NUM_BREAKPOINTS": BEGIN
         struct_is_ppoly = TRUE
         result = {DOMAIN_DIM: IMSL_LONG(sp.DOMAIN_DIM), $
           TARGET_DIM: IMSL_LONG(sp.TARGET_DIM), $
           ORDER: IMSL_LONG(sp.ORDER), $
           NUM_COEF: IMSL_LONG(sp.NUM_COEF), $
           NUM_BREAKPOINTS: IMSL_LONG(sp.NUM_BREAKPOINTS), $
           BREAKPOINTS: FIX(sp.BREAKPOINTS, TYPE=type), $
           COEF: FIX(sp.COEF, TYPE=type)}
      END
      "NUM_KNOTS": BEGIN
         struct_is_ppoly = FALSE
         result = {DOMAIN_DIM: IMSL_LONG(sp.DOMAIN_DIM), $
           TARGET_DIM: IMSL_LONG(sp.TARGET_DIM), $
           ORDER: IMSL_LONG(sp.ORDER), $
           NUM_COEF: IMSL_LONG(sp.NUM_COEF), $
           NUM_KNOTS: IMSL_LONG(sp.NUM_KNOTS), $
           KNOTS: FIX(sp.KNOTS, TYPE=type), $
           COEF: FIX(sp.COEF, TYPE=type)}
      END
      ELSE: message, err_str
   END

   ;Check the size of sp.(KNOTS/BREAKPOINTS) and sp.(COEF)
   total_coefs = IMSL_1
   FOR i = 0, sp.domain_dim-1 DO total_coefs = total_coefs*(sp.num_coef)(i)
   exp_sizes = [total(sp.(4)), sp.target_dim*total_coefs]
   IF (N_ELEMENTS(sp.(5)) NE exp_sizes(0)) THEN message, err_str
   IF (N_ELEMENTS(sp.(6)) NE exp_sizes(1)) THEN message, err_str
   RETURN, result
END
