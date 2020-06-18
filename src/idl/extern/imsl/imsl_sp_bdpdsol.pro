; $Id: //depot/idl/IDL_70/idl_src/libs/imsl/imsl_6.0/lib/imsl_sp_bdpdsol.pro#1 $
;
; Copyright (c) 1970-2006, VISUAL NUMERICS Inc. All Rights Reserved.
; This software is confidential information which is proprietary to and a
; trade secret of Visual Numerics Inc.  Use, duplication or disclosure is
; subject to the terms of an appropriate license agreement.
;
FUNCTION imsl_sp_bdpdsol, b, $                         ;INPUT 1-D array: floating point
                   ncoda, $                       ;INPUT Scalar LONG
                   a, $                           ;INPUT 1-D array: floating point
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
   ; CASE 1: Two (2) positional arguments:
   ;  In this case the RHS is supplied through the positional
   ;  argument, and the factored system must be supplied through
   ;  the keyword FACTOR. In this case the following checks 
   ;  are performed.
   ;  - b must be a 1-D array. n = N_ELEMENTS(b)
   ;  - FACTOR must be supplied.
   ;  - FACTOR must be a 2-D array of size ((ncoda+1) x n).
   ;  - CONDITION can't be present.
   ;   
   ; CASE 2: Three (3) positional arguments:
   ;  In this case, the positional arguments the RHS, b, the number of 
   ;  upper codiagonals ncoda, AND the band matrix, A.
   ;  stored as a 1-D array.
   ;  - b must be a 1-D array. n = N_ELEMENTS(b)
   ;  - A must be an array of length (ncoda+1)*n
   ;
   ; GENERAL CHECKS:
   ; The following checks are made regardless of the number of positional
   ; arguments.
   ;  - ncoda must be GE 0 and LE N_ELEMENTS(B)
   ;
   nargs = n_params()
   IF ((nargs NE 2) AND (nargs NE 3)) THEN $
     message, 'Incorrect number of arguments.'     
   type_a = IMSL_0
   size_b = IMSL_SIZE(b)
   IF (size_b(0) NE 1) THEN  $
     message, 'B must be a 1-D array.'
   n_rows = IMSL_LONG(size_b(3))
   factored = IMSL_0
   ncoda_cvt = IMSL_LONG(ncoda(0))
   
   IF ((ncoda_cvt LT 0) OR (ncoda_cvt GT n_rows)) $
     THEN message, 'NCODA must be greater than or equal to zero ' + $
                   'and less than N_ELEMENTS(B).'
   IF (nargs EQ 2) THEN BEGIN
      IF (KEYWORD_SET(FACTOR) NE 1) THEN BEGIN
         message, 'The keyword FACTOR is  required ' + $
           'when only two positional arguments are passed.'
      END ELSE BEGIN
         size_factor   = IMSL_SIZE(factor)
         type_factor = size_factor(N_ELEMENTS(size_factor)-2)
         ldfac = ncoda_cvt+1
         IF (size_factor(0) NE 2) THEN message, 'FACTOR must be a 2-D array.'
         IF ((size_factor(1) NE ldfac) OR (size_factor(2) NE n_rows)) $
           THEN message, 'FACTOR is not the correct size.'
      END
      IF (ARG_PRESENT(CONDITION)) THEN $
        message, 'CONDITION is not valid for this usage of SP_BDPDSOL.'
      a = 0.0                   ; Just a place holder.
      factored = IMSL_1
      solve_only = IMSL_1
   END ELSE BEGIN ;  Three positional arguments.
      size_a = IMSL_SIZE(a)
      type_a = size_a(N_ELEMENTS(size_a)-2)
      IF ((size_a(0) NE 1) OR (N_ELEMENTS(a) NE  (n_rows*(ncoda_cvt+1)))) $
          THEN message, 'A must be a 1-D array of length  (N_ELEMENTS(B))*(NCODA+1).'
   END
   ;
   ; Decide on what precision to use.
   ;
   type  =  TYP_FLOAT
   IF (size_b(N_ELEMENTS(size_b)-2) EQ  TYP_DOUBLE) THEN type = TYP_DOUBLE
   IF (size_b(N_ELEMENTS(size_b)-2) EQ  TYP_COMPLEX) THEN type = TYP_COMPLEX
   ; If FACTOR is supplied, then check its precision.
   IF (nargs EQ 2) THEN BEGIN
      IF (size_factor(N_ELEMENTS(size_factor)-2) EQ  TYP_DOUBLE) THEN type = (type > TYP_DOUBLE)
      IF (size_factor(N_ELEMENTS(size_factor)-2) EQ  TYP_COMPLEX) THEN type = (type > TYP_COMPLEX)
   END ELSE BEGIN
      IF (size_a(N_ELEMENTS(size_a)-2) EQ  TYP_DOUBLE) THEN type = (type > TYP_DOUBLE)
      IF (size_a(N_ELEMENTS(size_a)-2) EQ  TYP_COMPLEX) THEN type = (type > TYP_COMPLEX)
   END
   IF (KEYWORD_SET(double) EQ true) THEN begin
      IF (type EQ TYP_FLOAT) THEN type = TYP_DOUBLE
      IF (type EQ TYP_COMPLEX) THEN type = TYP_DCMPLX
   END
   ;
   ; Setup the parameters for the call to the system function.
   ; Note, there is some code here for future support of TYP_COMPLEX.
   ;
   cmplx_scale = IMSL_1
   IF ((type EQ TYP_DCMPLX) OR (type EQ  TYP_COMPLEX)) THEN cmplx_scale = IMSL_2
   IF ((type EQ TYP_DOUBLE) OR (type EQ TYP_DCMPLX)) THEN BEGIN
      ; Result
      ;
      result = dblarr(cmplx_scale*n_rows)
      ; 
      ; Output
      IF (ARG_PRESENT(condition)) THEN condition_spc = double(0.0)
   END ELSE BEGIN
      ; Result
      ;
      result = fltarr(cmplx_scale*n_rows)
      ; 
      ; Output
      IF (ARG_PRESENT(condition)) THEN condition_spc = float(0.0)
   END
   b_cvt = imsl_cvt_arr(b, type)
   IF (nargs EQ 3) THEN $
      a_cvt = imsl_cvt_arr(a, type)
   IF (nargs EQ 2) THEN $
     factor_cvt = imsl_cvt_arr(factor, type)
   ;
   ; Call the system function.
   ;
   err_status = 0L
   MATHSTAT_203,  type, err_status, b_cvt, a_cvt, $
                   n_rows, ncoda_cvt,  $
                   solve_only, $
                   factor_cvt, $
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
