; $Id: //depot/idl/IDL_70/idl_src/libs/imsl/imsl_6.0/lib/imsl_chnndsol.pro#1 $      
;
; Copyright (c) 1970-2006, VISUAL NUMERICS Inc. All Rights Reserved.
; This software is confidential information which is proprietary to and a
; trade secret of Visual Numerics Inc.  Use, duplication or disclosure is
; subject to the terms of an appropriate license agreement.
;
FUNCTION imsl_chnndsol, b, $                  ;INPUT 1-D array: floating point
                a, $                       ;INPUT 2-D array: floating point 
                double=double, $           ;INPUT Scalar ON/OFF flag
                tolerance=tolerance, $     ;INPUT Scalar floating point
                inverse=inverse, $         ;OUTPUT 2-D array: floating point
                factor=factor              ;INPUT 2-D array: floating point

@imsl_init.pro
   ON_ERROR, on_err_action
   ;
   ; Error checking.
   ;  The two basic uses of this function are divided between the
   ;  cases when there are either 1 or two position arguments.
   ;  Case 1: Only one positional argument.
   ;          In this case the RHS is supplied through the positional
   ;          argument, and the factored system must be supplied through
   ;          the keyword FACTOR. In this case the following checks 
   ;          are performed.
   ;          - B must be a 1-D array. n = N_ELEMENTS(b)
   ;          - FACTOR must be supplied.
   ;          - FACTOR must be a 2-D square array of order n.
   ;          
   ;  Case 2: Two positional arguments.
   ;          In this case, both the RHS and the original system have
   ;          have been supplied. In this case the following checks are 
   ;          performed.
   ;          - b must be a 1-D array. n = N_ELEMENTS(b)
   ;          - A must be a 2-D square array of order n.
   ;          - FACTOR can't be present.
   ;    
   ;    
   ;          
   nargs = n_params()
   IF ((nargs NE 1) AND (nargs NE 2)) THEN $
     message, 'Incorrect number of arguments.'
   size_b = IMSL_SIZE(b)
   IF (size_b(0) NE 1) THEN message, 'B must be a 1-D array.'
   n = IMSL_N_ELEMENTS(b)
   IF (nargs EQ 1) THEN BEGIN
      IF (NOT KEYWORD_SET(factor)) THEN $
        message, 'FACTOR is required for this usage of CHNNDSOL.'
      size_factor = IMSL_SIZE(factor)
      IF (size_factor(0) NE 2) THEN message, 'FACTOR must be a 2-D square array.'
      IF ((size_factor(1) NE n) OR (size_factor(2) NE n)) $
        THEN message, 'FACTOR is not the correct size.'
   END ELSE BEGIN
      size_a = IMSL_SIZE(a)
      IF (size_a(0) NE 2) THEN message, 'A must be a 2-D square array.'
      IF ((size_a(1) NE n) OR (size_a(2) NE n)) $
        THEN message, 'A is not the correct size.'
      IF (KEYWORD_SET(FACTOR)) THEN $
        message, 'FACTOR is not valid for this usage of CHNNDFAC.'
   END
   ;
   ; Decide on what precision to use.
   ;
   type = TYP_FLOAT
   type = TYP_FLOAT
   IF (size_b(N_ELEMENTS(size_b)-2) EQ  TYP_DOUBLE) THEN type = TYP_DOUBLE
   IF (nargs EQ 2) THEN $
     IF (size_a(N_ELEMENTS(size_a)-2) EQ  TYP_DOUBLE) THEN type = TYP_DOUBLE
   IF (KEYWORD_SET(double) EQ true) THEN type = TYP_DOUBLE
   ;
   ; Setup the parameters for the call to the system function.
   ;
   ; Floating point arguments and keywords
   IF (type EQ TYP_DOUBLE) THEN BEGIN
      result = dblarr(n)
      ; 
      ; Input
      b_cvt = double(b)
      IF (nargs EQ 2) THEN a_cvt = double(transpose(a))
      IF (KEYWORD_SET(factor)) THEN factor_cvt = double(transpose(factor))
      ; 
      ; Output
      IF (ARG_PRESENT(inverse)) THEN inverse_spc = dblarr(n, n)
      IF (KEYWORD_SET(tolerance)) THEN tolerance_cvt = double(tolerance(0))
   END ELSE BEGIN
      result = fltarr(n)
      ; 
      ; Input
      b_cvt = float(b)
      IF (nargs EQ 2) THEN a_cvt = float(transpose(a))
      IF (KEYWORD_SET(factor)) THEN factor_cvt = float(transpose(factor))
      ; 
      ; Output
      IF (ARG_PRESENT(inverse)) THEN inverse_spc = fltarr(n, n)
      IF (KEYWORD_SET(tolerance)) THEN tolerance_cvt = float(tolerance(0))
   END
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
                              tolerance_cvt, $
                              IMSL_2, $
                              result
   ;
   ; Now copy over all output keywords results.
   IF (arg_present(inverse)) THEN inverse = transpose(inverse_spc)
   ;
   ; return
   RETURN, result
END

                   
                   
                   

  
      

  
