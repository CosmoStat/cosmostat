; $Id: //depot/idl/IDL_70/idl_src/libs/imsl/imsl_6.0/lib/imsl_sp_bdpdfac.pro#1 $
;
; Copyright (c) 1970-2006, VISUAL NUMERICS Inc. All Rights Reserved.
; This software is confidential information which is proprietary to and a
; trade secret of Visual Numerics Inc.  Use, duplication or disclosure is
; subject to the terms of an appropriate license agreement.
;
FUNCTION imsl_Sp_bdpdfac, a, $                         ;INPUT 1-D array: floating point
                   n, $                           ;INPUT Scalar LONG
                   ncoda, $                       ;INPUT Scalar LONG
                   double=double, $               ;INPUT Scalar ON/OFF flag
                   condition=condition            ;OUTPUT Scalar floating point
@imsl_init.pro
   ON_ERROR, on_err_action
   ;
   ; Error checking.
   ; Checks on positional args:
   ;  - A must be an array of length (ncoda+1)*n
   ;  - n must be GE 0
   ;  - ncoda must be GE 0 and LT n
   ;
   nargs = n_params()
   IF (nargs NE 3) THEN $
     message, 'Incorrect number of arguments.'     
   type_a = IMSL_0

   n_rows = IMSL_LONG(n(0))
   ncoda_cvt = IMSL_LONG(ncoda(0))
   
   IF (n_rows LT 0) $
     THEN message, 'N must be greater than or equal to zero.'
   IF ((ncoda_cvt LT 0) OR (ncoda_cvt GE n_rows)) $
     THEN message, 'NCODA must be greater than or equal to zero ' + $
                   'and less than N_ELEMENTS(B).'
   size_a = IMSL_SIZE(a)
   type_a = size_a(N_ELEMENTS(size_a)-2)
   IF ((size_a(0) NE 1) OR (N_ELEMENTS(a) NE  (n_rows*(ncoda_cvt+1)))) $
     THEN message, 'A must be a 1-D array of length N*(NCODA+1).'
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
   ; Floating point arguments and keywords
   ldfac = (ncoda_cvt+1)
   cmplx_scale = IMSL_1
   IF ((type EQ TYP_DCMPLX) OR (type EQ  TYP_COMPLEX)) THEN cmplx_scale = IMSL_2
   IF ((type EQ TYP_DOUBLE) OR (type EQ TYP_DCMPLX)) THEN BEGIN
      ; 
      ; Output
      result = dblarr(cmplx_scale*n_rows, ldfac)
      IF (ARG_PRESENT(condition)) THEN condition_spc = double(0.0)
   END ELSE BEGIN
      ; 
      ; Output
      result = fltarr(cmplx_scale*n_rows, ldfac)
      IF (ARG_PRESENT(condition)) THEN condition_spc = float(0.0)
   END
   ;
   ; Convert A
   ;
   a_cvt = imsl_cvt_arr(a, type)
   ;
   ;
   ; Call the system function.
   ;
   err_status = 0L
   MATHSTAT_202,  type, err_status, a_cvt, $
                   n_rows, ncoda_cvt,  $
                   condition_spc, $
                   result
 
   ; Now copy over all output keywords results.
   ;
   IF (ARG_PRESENT(condition) EQ TRUE) THEN $
     condition = condition_spc
   ;
   ; Return.
   ;
   IF (cmplx_scale EQ 2) THEN BEGIN 
      RETURN, imsl_cvt_arr(result, /back)
   END ELSE BEGIN
      RETURN, transpose(result)
   END

END
