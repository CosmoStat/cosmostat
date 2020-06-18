; $Id: //depot/idl/IDL_70/idl_src/libs/imsl/imsl_6.0/lib/imsl_sp_bdfac.pro#1 $
;
; Copyright (c) 1970-2006, VISUAL NUMERICS Inc. All Rights Reserved.
; This software is confidential information which is proprietary to and a
; trade secret of Visual Numerics Inc.  Use, duplication or disclosure is
; subject to the terms of an appropriate license agreement.
;
PRO imsl_sp_bdfac,      nlca, $                        ;INPUT Scalar LONG
                   nuca, $                        ;INPUT Scalar LONG
                   n_rows, $                      ;INPUT Scalar LONG
                   a, $                           ;INPUT 1-D array: floating point
                   pivot,       $                 ;OUTPUT 1-D array: LONG
                   factor,        $                  ;OUTPUT 2-D array: floating point
                   blk_factor=blk_factor, $       ;INPUT Scalar LONG
                   double=double, $               ;INPUT Scalar ON/OFF flag
                   condition=condition            ;OUTPUT Scalar floating point

@imsl_init.pro
   ON_ERROR, on_err_action
   ;
   ; Error checking.
   ; Checks on positional args:
   ; There is one case for the number of positional argument
   ; passed to this function.
   ; CASE 1: Six (6) positional arguments:
   ;  - n_rows must be GE 1
   ;  - nlca must be GE 0
   ;  - nuca must be GE 0
   ;  - A must be an array of length (nlca+nuca+1)*n_rows
   ;
   nargs = n_params()
   IF (nargs NE 6) THEN $
     message, 'Incorrect number of arguments.' 
   type_a = IMSL_0
   n_rows_cvt = IMSL_LONG(n_rows(0))
   IF (n_rows_cvt LT 1) THEN message, 'N_ROWS must be greater than 0.'
   
   nlca_cvt = IMSL_LONG(nlca(0))
   IF ((nlca_cvt LT 0) OR (nlca_cvt GT n_rows_cvt)) $
     THEN message, 'NLCA must be greater than or equal to zero ' + $
                    'and less than N_ROWS.'
   nuca_cvt = IMSL_LONG(nuca(0))
   IF ((nuca_cvt LT 0) OR (nuca_cvt GT n_rows_cvt)) $
     THEN message, 'NUCA must be greater than or equal to zero ' + $
                   'and less than N_ROWS.'
      
   size_a = IMSL_SIZE(a)
   type_a = size_a(N_ELEMENTS(size_a)-2)
   IF ((size_a(0) NE 1) OR (N_ELEMENTS(a) NE  (n_rows_cvt*(nlca_cvt+nuca_cvt+1)))) $
     THEN message, 'A must be a 1-D array of length  (N_ELEMENTS(B))*(NLCA+NUCA+1).'
   
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
   ; Input LONG argument(s)
   IF (KEYWORD_SET(blk_factor) EQ TRUE) THEN $
     blk_factor_cvt = IMSL_LONG(blk_factor(0))
   ; Output LONG argument(s)
   pivot_spc = IMSL_LONARR(n_rows_cvt)
   ;
   ; Floating point arguments and keywords
   ldfac = (IMSL_2*nlca_cvt+nuca_cvt+1)
   cmplx_scale = IMSL_1
   IF ((type EQ TYP_DCMPLX) OR (type EQ  TYP_COMPLEX)) THEN cmplx_scale = IMSL_2
   IF ((type EQ TYP_DOUBLE) OR (type EQ TYP_DCMPLX)) THEN BEGIN
      ; 
      ; Output
      factor_spc = dblarr(cmplx_scale*n_rows_cvt, ldfac)
      IF (ARG_PRESENT(condition)) THEN condition_spc = double(0.0)
   END ELSE BEGIN
      ; 
      ; Output
      factor_spc = fltarr(cmplx_scale*n_rows_cvt, ldfac)
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
   MATHSTAT_201,LONG(type), err_status, a_cvt, n_rows_cvt, nlca_cvt, nuca_cvt,  $
                   factor_spc, $
                   pivot_spc, $
                   blk_factor_cvt, $ 
                   condition_spc

   ; Now copy over all output keywords results.
   ;
   IF (ARG_PRESENT(condition)) THEN condition  =  condition_spc
   pivot = pivot_spc
   IF (cmplx_scale EQ 2) THEN BEGIN 
      factor = imsl_cvt_arr(factor_spc, /back)
   END ELSE BEGIN
      factor = transpose(factor_spc)
   END
    IF (ARG_PRESENT(condition) EQ TRUE) THEN $
     condition = condition_spc
   ;
   ; Return.
   ;
END

                   
                   
                   

  
      

  
