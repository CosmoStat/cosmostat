; $Id: //depot/idl/IDL_70/idl_src/libs/imsl/imsl_6.0/lib/imsl_sp_bdsol.pro#1 $
;
; Copyright (c) 1970-2006, VISUAL NUMERICS Inc. All Rights Reserved.
; This software is confidential information which is proprietary to and a
; trade secret of Visual Numerics Inc.  Use, duplication or disclosure is
; subject to the terms of an appropriate license agreement.
;
FUNCTION imsl_sp_bdsol, b, $                           ;INPUT 1-D array: floating point
                   nlca, $                        ;INPUT Scalar LONG
                   nuca, $                        ;INPUT Scalar LONG
                   a, $                           ;INPUT 1-D array: floating point
                   transpose=transpose, $         ;INPUT Scalar LONG
                   blk_factor=blk_factor, $       ;INPUT Scalar LONG
                   pivot=pivot, $                 ;INPUT 1-D array: LONG
                   factor=factor, $               ;INPUT 2-D array: floating point
                   double=double, $               ;INPUT Scalar ON/OFF flag
                   condition=condition            ;OUTPUT Scalar floating point

@imsl_init.pro
   ON_ERROR, on_err_action
   ;
   ; Error checking.
   ; Checks on positional args:
   ; There are two different cases for the number of positional argument
   ; passed to this function.
   ; CASE 1: Three (3) positional arguments:
   ;  In this case the RHS is supplied through the positional
   ;  argument, and the factored system must be supplied through
   ;  the keyword FACTOR. In this case the following checks 
   ;  are performed.
   ;  - b must be a 1-D array. n = N_ELEMENTS(b)
   ;  - FACTOR and PIVOT must be supplied.
   ;  - FACTOR must be a 2-D array of size ((2*nlca+nuca+1) x n).
   ;  - PIVOT must be a 1-D array of length n.
   ;  - CONDITION can't be present.
   ;   
   ; CASE 2: Four (4) positional arguments:
   ;  In this case, the positional arguments the RHS, b, the number of 
   ;  codiagonals nlca & nuca, AND the band matrix, A.
   ;  stored as a 1-D array, and .
   ;  - b must be a 1-D array. n = N_ELEMENTS(b)
   ;  - A must be an array of length (nlca+nuca+1)*n_rows
   ;
   ; GENERAL CHECKS:
   ; The following checks are made regardless of the number of positional
   ; arguments.
   ;  - nlca must be GE 0 and LE n
   ;  - nuca must be GE 0 and LE n
   ;
   nargs = n_params()
   IF ((nargs NE 3) AND (nargs NE 4)) THEN $
     message, 'Incorrect number of arguments.'     
   type_a = IMSL_0
   size_b = IMSL_SIZE(b)
   IF (size_b(0) NE 1) THEN  $
     message, 'B must be a 1-D array.'
   n_rows = IMSL_LONG(size_b(3))
   factored = IMSL_0
   nlca_cvt = IMSL_LONG(nlca(0))
   IF ((nlca_cvt LT 0) OR (nlca_cvt GT n_rows)) $
     THEN message, 'NLCA must be greater than or equal to zero ' + $
                   'and less than N_ELEMENTS(B).'
   nuca_cvt = IMSL_LONG(nuca(0))
   IF ((nuca_cvt LT 0) OR (nuca_cvt GT n_rows)) $
     THEN message, 'NUCA must be greater than or equal to zero ' + $
                   'and less than N_ELEMENTS(B).'
   IF (nargs EQ 3) THEN BEGIN
      IF ((KEYWORD_SET(PIVOT) + KEYWORD_SET(FACTOR)) NE 2) THEN BEGIN
         message, 'The keywords PIVOT and FACTOR are required ' + $
           'when only three positional arguments are passed.'
      END ELSE BEGIN
         size_pivot     = IMSL_SIZE(pivot)
         size_factor   = IMSL_SIZE(factor)
         
         type_factor = size_factor(N_ELEMENTS(size_factor)-2)
         IF ((size_pivot(0) NE 1) OR (N_ELEMENTS(pivot) NE (n_rows))) THEN $
           message, 'PIVOT must be a 1-D array OF length N_ELEMENTS(B).'
         pivot_cvt = IMSL_LONG(pivot)
         
         ldfac = IMSL_2*nlca_cvt+nuca_cvt+1
         IF (size_factor(0) NE 2) THEN message, 'FACTOR must be a 2-D array.'
         IF ((size_factor(1) NE ldfac) OR (size_factor(2) NE n_rows)) $
           THEN message, 'FACTOR is not the correct size.'
      END
      IF (ARG_PRESENT(CONDITION)) THEN $
        message, 'CONDITION is not valid for this usage of SP_BDSOL.'
      a = 0.0                   ; Just a place holder.
      factored = IMSL_1
      solve_only = IMSL_1
   END ELSE BEGIN ;  Four positional arguments.
      
      size_a = IMSL_SIZE(a)
      type_a = size_a(N_ELEMENTS(size_a)-2)
      IF ((size_a(0) NE 1) OR (N_ELEMENTS(a) NE  (n_rows*(nlca_cvt+nuca_cvt+1)))) $
          THEN message, 'A must be a 1-D array of length  (N_ELEMENTS(B))*(NLCA+NUCA+1).'
   END
   ;
   ; Decide on what precision to use.
   ;
   type  =  TYP_FLOAT
   IF (size_b(N_ELEMENTS(size_b)-2) EQ  TYP_DOUBLE) THEN type = TYP_DOUBLE
   IF (size_b(N_ELEMENTS(size_b)-2) EQ  TYP_COMPLEX) THEN type = TYP_COMPLEX
   IF (size_b(N_ELEMENTS(size_b)-2) EQ TYP_DCMPLX) THEN type = TYP_DCMPLX
   ; If FACTOR is supplied, then check its precision.
   IF (nargs EQ 3) THEN BEGIN
      IF (size_factor(N_ELEMENTS(size_factor)-2) EQ  TYP_DOUBLE) THEN type = (type > TYP_DOUBLE)
      IF (size_factor(N_ELEMENTS(size_factor)-2) EQ  TYP_COMPLEX) THEN type = (type > TYP_COMPLEX)
      IF (size_factor(N_ELEMENTS(size_factor)-2) EQ TYP_DCMPLX) THEN type = (type > TYP_DCMPLX)
   END ELSE BEGIN
      IF (size_a(N_ELEMENTS(size_a)-2) EQ  TYP_DOUBLE) THEN type = (type > TYP_DOUBLE)
      IF (size_a(N_ELEMENTS(size_a)-2) EQ  TYP_COMPLEX) THEN type = (type > TYP_COMPLEX)
      IF (size_a(N_ELEMENTS(size_a)-2) EQ TYP_DCMPLX) THEN type = (type > TYPDCMPLX)
   END
   
   IF (KEYWORD_SET(double) EQ true) THEN begin
      IF (type EQ TYP_FLOAT) THEN type = TYP_DOUBLE
      IF (type EQ TYP_COMPLEX) THEN type = TYP_DCMPLX
   END
   ;
   ; Setup the parameters for the call to the system function.
   ; Note, there is some code here for future support of TYP_COMPLEX.
   ;
   ; Input LONG argument(s)
   IF (KEYWORD_SET(transpose) EQ TRUE) THEN $
     transpose_cvt = IMSL_1
   IF (KEYWORD_SET(blk_factor) EQ TRUE) THEN $
     blk_factor_cvt = IMSL_LONG(blk_factor(0))
   IF (KEYWORD_SET(pivot) EQ TRUE) THEN $ 
     pivot_cvt = IMSL_LONG(pivot)
   cmplx_scale = IMSL_1
   IF ((type EQ TYP_DCMPLX) OR (type EQ  TYP_COMPLEX)) THEN cmplx_scale = IMSL_2

   ;
   ; Floating point arguments and keywords
   ;
   IF ((type EQ TYP_DOUBLE) OR (type EQ TYP_DCMPLX)) THEN BEGIN
      result = dblarr(cmplx_scale*n_rows)
      ; 
      ; Output
      IF (ARG_PRESENT(condition)) THEN condition_spc = double(0.0)
   END ELSE BEGIN
      result = fltarr(cmplx_scale*n_rows)
      ; 
      ; Output
      IF (ARG_PRESENT(condition)) THEN condition_spc = float(0.0)
   END
   IF (nargs EQ 4) THEN $
      a_cvt = imsl_cvt_arr(a, type)
   b_cvt = imsl_cvt_arr(b, type)
   IF (nargs EQ 3) THEN $
     factor_cvt = imsl_cvt_arr(factor, type)
   ;
   ; Call the system function.
   ;
   err_status = 0L
   MATHSTAT_204, LONG(type), err_status, b_cvt, a_cvt, $
                   n_rows, nlca_cvt, nuca_cvt,  $
                   transpose_cvt, $
                   solve_only, $
                   factor_cvt, $
                   pivot_cvt, $
                   blk_factor_cvt, $ 
                   condition_spc, $
                   result

   ; Now copy over all output keywords results.
   ;
   IF (ARG_PRESENT(condition) EQ TRUE) THEN $
     condition = condition_spc
   ;
   ; Return.
   ;
   IF (cmplx_scale EQ 2) $
     THEN RETURN, imsl_cvt_arr(result, /back) $
     ELSE RETURN, result
END

                   
                   
                   

  
      

  
