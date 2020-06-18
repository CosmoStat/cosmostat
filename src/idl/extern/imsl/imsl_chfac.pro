; $Id: //depot/idl/IDL_70/idl_src/libs/imsl/imsl_6.0/lib/imsl_chfac.pro#1 $      
;
; Copyright (c) 1970-2006, VISUAL NUMERICS Inc. All Rights Reserved.
; This software is confidential information which is proprietary to and a
; trade secret of Visual Numerics Inc.  Use, duplication or disclosure is
; subject to the terms of an appropriate license agreement.
;
PRO      imsl_chfac, a, $                     ;INPUT 2-D array: floating point 
                fac, $                     ;OUTPUT 1-D array: LONG
                complex=complex, $         ;INPUT Scalar ON/OFF flag
                double=double, $           ;INPUT Scalar ON/OFF flag
                inverse=inverse, $         ;OUTPUT 2-D array: floating point
                condition=condition        ;OUTPUT Scalar floating point

@imsl_init.pro
   ON_ERROR, on_err_action
   ;
   ; Error checking:
   ;  The two basic uses of this function are divided between the
   ;  cases when there are either 1 two, or three position arguments.
   ;  Case 1: Only one positional argument.
   ;          In this case the factor will not be returned
   ;          through a positional argument.  
   ;          In this case the following checks are performed.
   ;          - A must be a 2-D square array.
   ;          
   ;  Case 2: Two positional arguments.
   ;          In this case the factor will be returned
   ;          through the a positional argument.  
   ;          In this case the following checks are performed.
   ;          - A must be a 2-D square array.
   ;          
   ;    
   ;   Later on in the code we check that the keyword INVERSE is
   ;   is not used if the returned data type is TYP_COMPLEX.
   ;          
   ;          
   nargs = n_params()
   IF ((nargs LT 1) OR (nargs GT 2)) THEN $
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
   IF (size_a(N_ELEMENTS(size_a)-2) EQ  TYP_COMPLEX) THEN type = TYP_COMPLEX
   IF (size_a(N_ELEMENTS(size_a)-2) EQ TYP_DCMPLX) THEN type = TYP_DCMPLX
   IF (KEYWORD_SET(double) EQ true) THEN begin
      IF (type EQ TYP_FLOAT) THEN type = TYP_DOUBLE
      IF (type EQ TYP_COMPLEX) THEN type = TYP_DCMPLX
   END
   IF (KEYWORD_SET(complex) EQ true) THEN begin
      IF (type EQ TYP_FLOAT) THEN type = TYP_COMPLEX
      IF (type EQ TYP_DOUBLE) THEN type = TYP_DCMPLX
   END
   IF (KEYWORD_SET(inverse) AND $
       ((type EQ TYP_COMPLEX) OR (type EQ TYP_DCMPLX))) THEN $
     message, 'FACTOR is not valid for this usage of CHFAC.'
   ;
   ; Setup the parameters for the call to the system function.
   ;
   ; Floating point arguments and keywords
   cmplx_scale = IMSL_1
   IF ((type EQ TYP_DCMPLX) OR (type EQ  TYP_COMPLEX)) THEN cmplx_scale = IMSL_2
   IF ((type EQ TYP_DOUBLE) OR (type EQ TYP_DCMPLX)) THEN BEGIN
      ; 
      ; Output
      fac_spc = dblarr(cmplx_scale*n, n)
      IF (ARG_PRESENT(inverse)) THEN inverse_spc = dblarr(cmplx_scale*n, n)
      IF (ARG_PRESENT(condition)) THEN condition_spc = double(0.0)
   END ELSE BEGIN
      ; 
      ; Output
      fac_spc = fltarr(cmplx_scale*n, n)
      IF (ARG_PRESENT(inverse)) THEN inverse_spc = fltarr(cmplx_scale*n, n)
      IF (ARG_PRESENT(condition)) THEN condition_spc = float(0.0)
   END
     
   ;
   ; Convert A
   ;
   a_cvt = imsl_cvt_arr(a, type)
   ;
   ; Call the system function.
   ;
   err_status = 0L
   MATHSTAT_157, type, err_status, a_cvt, fac_spc, n, $
                              pivot_spc, $      ;NOT USED
                              transpose_cvt, $  ;NOT USED
                              inverse_spc, $
                              condition_spc, $
                              tolerance_cvt, $  ;NOT USED
                              IMSL_1
   ;
   ; Now copy over all output keywords results.
   IF (arg_present(condition)) THEN condition = condition_spc
   IF (cmplx_scale EQ 2) THEN BEGIN 
      IF (arg_present(inverse)) THEN inverse = imsl_cvt_arr(inverse_spc, /back)
      fac = imsl_cvt_arr(fac_spc, /back)
   END ELSE BEGIN
      IF (arg_present(inverse)) THEN inverse = transpose(inverse_spc)
      fac = transpose(fac_spc)
   END
   ;
   ; End of procedure.
END
