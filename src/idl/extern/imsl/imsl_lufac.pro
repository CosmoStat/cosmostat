; $Id: //depot/idl/IDL_70/idl_src/libs/imsl/imsl_6.0/lib/imsl_lufac.pro#1 $      
;
; Copyright (c) 1970-2006, VISUAL NUMERICS Inc. All Rights Reserved.
; This software is confidential information which is proprietary to and a
; trade secret of Visual Numerics Inc.  Use, duplication or disclosure is
; subject to the terms of an appropriate license agreement.
;
PRO      imsl_lufac, a, $                     ;INPUT 2-D array: floating point 
                pivot, $                   ;OUTPUT 1-D array: LONG
                fac, $                     ;OUTPUT 1-D array: LONG
                transpose=transpose, $     ;INPUT Scalar ON/OFF flag
                complex=complex, $         ;INPUT Scalar ON/OFF flag
                double=double, $           ;INPUT Scalar ON/OFF flag
                inverse=inverse, $         ;OUTPUT 2-D array: floating point
                condition=condition, $     ;OUTPUT Scalar floating point
                l=l, $                     ;OUTPUT 2-D array: floating point
                pa=pa, $                   ;OUTPUT 2-D array: floating point
                u=u                        ;OUTPUT 2-D array: floating point
@imsl_init.pro
   ON_ERROR, on_err_action
   ;
   ; Error checking:
   ;  The two basic uses of this function are divided between the
   ;  cases when there are either one, two, or three position arguments.
   ;  Case 1: Only one positional argument.
   ;          In this case the pivot and factor will not be returned
   ;          through positional arguments.  
   ;          In this case the following checks are performed.
   ;          - A must be a 2-D square array.
   ;          
   ;  Case 2: Two positional arguments.
   ;          In this case the pivot will be returned
   ;          through second positional argument.  
   ;          In this case the following checks are performed.
   ;          - A must be a 2-D square array.
   ;
   ;  Case 3: Three positional arguments.
   ;          In this case the pivot and factor will be returned
   ;          through the positional arguments.  
   ;          In this case the following checks are performed.
   ;          - A must be a 2-D square array.
   ;          
   ;          
   nargs = n_params()
   IF ((nargs LT 1) OR (nargs GT 3)) THEN $
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
   ;
   ; Setup the parameters for the call to the system function.
   ;
   ; Input LONG keyword(s)
   IF (KEYWORD_SET(transpose) EQ TRUE) THEN transpose_cvt = IMSL_1
   pivot_spc = IMSL_LONARR(n)
   ; Floating point arguments and keywords
   cmplx_scale = IMSL_1
   IF ((type EQ TYP_DCMPLX) OR (type EQ  TYP_COMPLEX)) THEN cmplx_scale = IMSL_2
   IF ((type EQ TYP_DOUBLE) OR (type EQ  TYP_DCMPLX)) THEN BEGIN
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
                              pivot_spc, $ 
                              transpose_cvt, $
                              inverse_spc, $
                              condition_spc, $
                              tolerance_cvt, $ ;NOT USED
                              IMSL_0
   ;
   ; Now copy over all output keywords results.
   IF (arg_present(condition)) THEN condition = condition_spc
   pivot = pivot_spc
   IF (cmplx_scale EQ 2) THEN BEGIN 
      IF (arg_present(inverse)) THEN inverse = imsl_cvt_arr(inverse_spc, /back)
      fac = imsl_cvt_arr(fac_spc, /back)
   END ELSE BEGIN
      IF (arg_present(inverse)) THEN inverse = transpose(inverse_spc)
      fac = transpose(fac_spc)
   END
   ;
   ; The keywords L, U, and PA are produced from the
   ; the contents of FAC, PIVOT and A.
   ; Compute L if needed
   IF (arg_present(l)) THEN BEGIN
      l = make_array(size = IMSL_SIZE(a_cvt))
      ; Copy the lower part of FAC into L, and put 1's on diagonal.
      FOR i = 1, n-1 DO l(i,0:i-1) = fac(i, 0:i-1)
      l(lindgen(n), lindgen(n)) = 1.0

      ;Unscramble L
      FOR k = 1, n-1 DO $
        IF (k NE (pivot(k)-1)) THEN BEGIN
         t = l(pivot(k)-1, 0:(k-1))
         l(pivot(k)-1, 0:(k-1)) = l(k, 0:(k-1))
         l(k, 0:(k-1)) = t
      END
      
      ; Correct signs of L
      FOR i = 1, n-1 DO l(i, 0:i-1) = -1.*l(i, 0:i-1)
   END
   ; 
   ; Compute U if needed
   IF (arg_present(u)) THEN BEGIN
      u = make_array(size = IMSL_SIZE(a_cvt))
      ; Copy the upper part of FAC into U.
      FOR i = 0, n-1 DO u(i, i:*) = fac(i, i:*)
   END
   ; 
   ; Compute PA if needed
   IF (arg_present(pa)) THEN BEGIN
      pa = make_array(size = IMSL_SIZE(a_cvt))
      pa = a
      FOR j = 0, n-1 DO BEGIN
         k = pivot(j) - 1
         IF (k NE j) THEN BEGIN
            t = pa(j, *)
            pa(j, *) = pa(k, *)
            pa(k, *) = t
         END
      END
   END
   ;
   ; End of procedure.
END
