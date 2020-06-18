; $Id: //depot/idl/IDL_70/idl_src/libs/imsl/imsl_6.0/lib/imsl_chsol.pro#1 $      
;
; Copyright (c) 1970-2006, VISUAL NUMERICS Inc. All Rights Reserved.
; This software is confidential information which is proprietary to and a
; trade secret of Visual Numerics Inc.  Use, duplication or disclosure is
; subject to the terms of an appropriate license agreement.
;
FUNCTION imsl_chsol, b, $                     ;INPUT 1-D array: floating point
                a, $                       ;INPUT 2-D array: floating point 
                complex=complex, $         ;INPUT Scalar ON/OFF flag
                double=double, $           ;INPUT Scalar ON/OFF flag
                inverse=inverse, $         ;OUTPUT 2-D array: floating point
                condition=condition, $     ;OUTPUT Scalar floating point
                factor=factor              ;OUTPUT 2-D array: floating point

@imsl_init.pro
   ON_ERROR, on_err_action
   ;
   ; Error checking.
   ;  The two basic uses of this function are divided between the
   ;  cases when there are either one or two position arguments.
   ;  Case 1: Only one positional argument.
   ;          In this case the RHS is supplied through the positional
   ;          argument, and the factored system must be supplied through
   ;          the keyword FACTOR. In this case the following checks 
   ;          are performed.
   ;          - b must be a 1-D array. n = N_ELEMENTS(b)
   ;          - FACTOR must be supplied.
   ;          - FACTOR must be a 2-D square array of order n.
   ;          - CONDITION can't be present.
   ;          
   ;  Case 2: Two positional arguments.
   ;          In this case, both the RHS and the original system have
   ;          have been supplied. In this case the following checks are 
   ;          performed.
   ;          - B must be a 1-D array. n = N_ELEMENTS(b)
   ;          - A must be a 2-D square array of order n.
   ;          - FACTOR can't be present.
   ;
   ;  Later on in the code we check that the keyword INVERSE is
   ;  is not used if the returned data type is TYP_COMPLEX.
   ;          
   nargs = n_params()
   IF ((nargs NE 1) AND (nargs NE 2)) THEN $
     message, 'Incorrect number of arguments.'
   size_b = IMSL_SIZE(b)
   IF (size_b(0) NE 1) THEN message, 'B must be a 1-D array.'
   n = IMSL_N_ELEMENTS(b)
   IF (nargs EQ 1) THEN BEGIN
      IF (NOT KEYWORD_SET(factor)) THEN $
        message, 'FACTOR is required for this usage of CHSOL.'
      size_factor = IMSL_SIZE(factor)
      IF (size_factor(0) NE 2) THEN message, 'FACTOR must be a 2-D square array.'
      IF ((size_factor(1) NE n) OR (size_factor(2) NE n)) $
        THEN message, 'FACTOR is not the correct size.'
      IF (ARG_PRESENT(CONDITION)) THEN $
        message, 'CONDITION is not valid for this usage of CHSOL.'
   END ELSE BEGIN
      size_a = IMSL_SIZE(a)
      IF (size_a(0) NE 2) THEN message, 'A must be a 2-D square array.'
      IF ((size_a(1) NE n) OR (size_a(2) NE n)) $
        THEN message, 'A is not the correct size.'
      IF (KEYWORD_SET(FACTOR)) THEN $
        message, 'FACTOR is not valid for this usage of CHSOL.'
   END
   ;
   ; Decide on what precision to use.
   ;
   type = TYP_FLOAT
   IF (size_b(N_ELEMENTS(size_b)-2) EQ  TYP_DOUBLE) THEN type = TYP_DOUBLE
   IF (size_b(N_ELEMENTS(size_b)-2) EQ  TYP_COMPLEX) THEN type = TYP_COMPLEX
   IF (size_b(N_ELEMENTS(size_b)-2) EQ TYP_DCMPLX) THEN type = TYP_DCMPLX
   IF (nargs EQ 2) THEN BEGIN
      IF (size_a(N_ELEMENTS(size_a)-2) EQ  TYP_DOUBLE) THEN type = TYP_DOUBLE
      IF (size_a(N_ELEMENTS(size_a)-2) EQ  TYP_COMPLEX) THEN type = TYP_COMPLEX
      IF (size_a(N_ELEMENTS(size_a)-2) EQ TYP_DCMPLX) THEN type = TYP_DCMPLX
   END
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
     message, 'FACTOR is not valid for this usage of CHSOL.'
   ;
   ; Setup the parameters for the call to the system function.
   ;
   ; Input LONG keyword(s)
   ; Floating point arguments and keywords
   cmplx_scale = IMSL_1
   IF ((type EQ TYP_DCMPLX) OR (type EQ  TYP_COMPLEX)) THEN cmplx_scale = IMSL_2
   IF ((type EQ TYP_DOUBLE) OR (type EQ TYP_DCMPLX)) THEN BEGIN
      result = dblarr(cmplx_scale*n)
      ; 
      ; Output
      IF (ARG_PRESENT(inverse)) THEN inverse_spc = dblarr(cmplx_scale*n, n)
      IF (ARG_PRESENT(condition)) THEN condition_spc = double(0.0)
   END ELSE BEGIN
      result = fltarr(cmplx_scale*n)
      ; 
      ; Output
      IF (ARG_PRESENT(inverse)) THEN inverse_spc = fltarr(cmplx_scale*n, n)
      IF (ARG_PRESENT(condition)) THEN condition_spc = float(0.0)
   END
     
   ;
   ; Convert B, A and FACTOR.
   ;
   b_cvt = imsl_cvt_arr(b, type)
   IF (nargs EQ 2) THEN a_cvt = imsl_cvt_arr(a, type)
   IF (KEYWORD_SET(factor)) THEN factor_cvt = imsl_cvt_arr(factor, type)
   ;
   ; Call the system function.
   ;
   err_status = 0L
   MATHSTAT_158, type, err_status, b_cvt, a_cvt, n, $
                              transpose, $ ;NOT USED
                              pivot, $     ;NOT USED
                              factor_cvt, $
                              inverse_spc, $
                              condition_spc, $
                              tolerance_cvt, $ ;NOT USED
                              IMSL_1, $
                              result
   ;
   ; Now copy over all output keywords results.
   IF (arg_present(condition)) THEN condition = condition_spc
   IF (arg_present(inverse)) THEN $
     IF (cmplx_scale EQ 2) THEN inverse = imsl_cvt_arr(inverse_spc, /back) $
     ELSE inverse = transpose(inverse_spc)
   ;
   ; return
   IF (cmplx_scale EQ 2) $
     THEN RETURN, imsl_cvt_arr(result, /back) $
     ELSE RETURN, result
END

                   
                   
                   

  
      

  
