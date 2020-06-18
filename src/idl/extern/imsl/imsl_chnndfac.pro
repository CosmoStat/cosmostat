; $Id: //depot/idl/IDL_70/idl_src/libs/imsl/imsl_6.0/lib/imsl_chnndfac.pro#1 $      
;
; Copyright (c) 1970-2006, VISUAL NUMERICS Inc. All Rights Reserved.
; This software is confidential information which is proprietary to and a
; trade secret of Visual Numerics Inc.  Use, duplication or disclosure is
; subject to the terms of an appropriate license agreement.
;
PRO      imsl_chnndfac, a, $                  ;INPUT 2-D array: floating point 
                fac, $                     ;OUTPUT 1-D array: LONG
                double=double, $           ;INPUT Scalar ON/OFF flag
                inverse=inverse, $         ;OUTPUT 2-D array: floating point
                tolerance=tolerance        ;INPUT Scalar floating point

@imsl_init.pro
   ON_ERROR, on_err_action
   ;
   ; Error checking:
   ;  - A must be a 2-D square array.
   ;          
   nargs = n_params()
   IF (nargs NE 2)  THEN $
     message, 'Incorrect number of arguments.'
   size_a = IMSL_SIZE(a)
   IF (size_a(0) NE 2) THEN message, 'A must be a 2-D square array.'
   IF (size_a(1) NE size_a(2)) $
     THEN message, 'A is not the correct size.'
   n = IMSL_LONG(size_a(1))
   ;
   ; Decide on what precision to use.
   ;
   type = TYP_FLOAT
   IF (size_a(N_ELEMENTS(size_a)-2) EQ  TYP_DOUBLE) THEN type = TYP_DOUBLE
   IF (KEYWORD_SET(double) EQ true) THEN type = TYP_DOUBLE
   ;
   ; Setup the parameters for the call to the system function.
   ;
   ; Floating point arguments and keywords
   IF (type EQ TYP_DOUBLE) THEN BEGIN
      ; 
      ; Input
      a_cvt = double(a)
      ; 
      ; Output
      fac_spc = dblarr(n, n)
      IF (ARG_PRESENT(inverse)) THEN inverse_spc = dblarr(n, n)
      IF (KEYWORD_SET(tolerance)) THEN tolerance_cvt = double(tolerance(0))
   END ELSE BEGIN
      ; 
      ; Input
      a_cvt = float(a)
      ; 
      ; Output
      fac_spc = fltarr(n, n)
      IF (ARG_PRESENT(inverse)) THEN inverse_spc = fltarr(n, n)
      IF (KEYWORD_SET(tolerance)) THEN tolerance_cvt = float(tolerance(0))
   END
   ;
   ; Call the system function.
   ;
   err_status = 0L
   MATHSTAT_157, type, err_status, a_cvt, fac_spc, n, $
                              pivot_spc, $      ;NOT USED
                              transpose_cvt, $  ;NOT USED
                              inverse_spc, $
                              condition_spc, $  ;NOT USED
                              tolerance_cvt, $
                              IMSL_2
   ;
   ; Now copy over all output keywords results.
   IF (arg_present(inverse)) THEN inverse = transpose(inverse_spc)
   fac = transpose(fac_spc)
   ;
   ; End of procedure.
END
