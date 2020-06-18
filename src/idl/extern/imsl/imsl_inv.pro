; $Id: //depot/idl/IDL_70/idl_src/libs/imsl/imsl_6.0/lib/imsl_inv.pro#1 $   
;
; Copyright (c) 1970-2006, VISUAL NUMERICS Inc. All Rights Reserved.
; This software is confidential information which is proprietary to and a
; trade secret of Visual Numerics Inc.  Use, duplication or disclosure is
; subject to the terms of an appropriate license agreement.
;
   
FUNCTION imsl_inv, a, $                       ;INPUT 2-D array: floating point 
                double=double              ;INPUT Scalar ON/OFF flag

@imsl_init.pro
   ON_ERROR, on_err_action
   ;
   ; Error checking:
   ;    - A must be 2-D.
   ;    - A must be square

   nargs = n_params()
   IF (nargs NE 1) THEN message, 'Incorrect number of arguments.'
   size_a = IMSL_SIZE(a)
   IF (size_a(0) NE 2) THEN message, 'A must be a 2-D square array.'
   IF (size_a(1) NE size_a(2)) THEN message, 'A must be a 2-D square array.' $
     ELSE n = IMSL_LONG(size_a(1))
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
   ;
   ; Setup the parameters for the call to the system function.
   ;
   ; Floating point arguments and keywords
   cmplx_scale = IMSL_1
   IF ((type EQ TYP_DCMPLX) OR (type EQ  TYP_COMPLEX)) THEN cmplx_scale = IMSL_2
   IF ((type EQ TYP_FLOAT) OR (type EQ  TYP_COMPLEX)) THEN $
     result = fltarr(cmplx_scale*n, n) $
   ELSE $
     result = dblarr(cmplx_scale*n, n)
   ;
   ; Convert A.
   ;
   a_cvt = imsl_cvt_arr(a, type)
   ;
   ; Call the system function.
   ;
   err_status = 0L
   MATHSTAT_152, type, err_status, a_cvt, n, $ 
                              result
   ;
   ; Now copy over all output keywords results.
   ;
   ; return
   IF (cmplx_scale EQ 2) $
     THEN RETURN, imsl_cvt_arr(result, /back) $
     ELSE RETURN, transpose(result)
END

                   
                   
                   

  
      

  
