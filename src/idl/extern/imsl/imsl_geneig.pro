; $Id: //depot/idl/IDL_70/idl_src/libs/imsl/imsl_6.0/lib/imsl_geneig.pro#1 $
;
; Copyright (c) 1970-2006, VISUAL NUMERICS Inc. All Rights Reserved.
; This software is confidential information which is proprietary to and a
; trade secret of Visual Numerics Inc.  Use, duplication or disclosure is
; subject to the terms of an appropriate license agreement.
;
PRO imsl_geneig, a, $                           ;INPUT 2-D array: floating point 
                 b, $                       ;INPUT 2-D array: floating point 
                 alpha, $                   ;OUTPUT 1-D array: floating point 
                 beta, $                    ;OUTPUT 1-D array: floating point 
                 DOUBLE = DOUBLE,   $       ;INPUT Scalar ON/OFF flag
                 vectors=vectors            ;OUTPUT 2-D array: floating point 

@imsl_init.pro
   ON_ERROR, on_err_action
   ;
   ; Error checking:
   ;    - A must be 2-D.
   ;    - A must be square, set n based on A.
   ;    - B must be 2-D of size (n, n)

   ;
   nargs  =  N_PARAMS()
   IF (nargs NE 4) THEN message, 'Incorrect number of arguments.'
   size_a = IMSL_SIZE(a)
   IF (size_a(0) NE 2) THEN message, 'A must be a 2-D square array.'
   IF (size_a(1) NE size_a(2)) THEN message, 'A must be a 2-D square array.' $
     ELSE n = IMSL_LONG(size_a(1))
   size_b = IMSL_SIZE(b)
   IF (size_b(0) NE 2) THEN message, 'B must be a 2-D square array.'
   IF ((size_b(1) NE n) OR (size_b(2) NE n)) THEN message, 'B must be the same size as A.'

   ;
   ; Decide on what precision to use.
   ;
   type = TYP_FLOAT
   IF (size_a(N_ELEMENTS(size_a)-2) EQ  TYP_DOUBLE) THEN type = TYP_DOUBLE
   IF (size_a(N_ELEMENTS(size_a)-2) EQ  TYP_COMPLEX) THEN type = TYP_COMPLEX
   IF (size_a(N_ELEMENTS(size_a)-2) EQ TYP_DCMPLX) THEN type = TYP_DCMPLX
   IF (size_b(N_ELEMENTS(size_b)-2) EQ  TYP_DOUBLE) THEN type = TYP_DOUBLE
   IF (size_b(N_ELEMENTS(size_b)-2) EQ  TYP_COMPLEX) THEN type = TYP_COMPLEX
   IF (size_b(N_ELEMENTS(size_b)-2) EQ TYP_DCMPLX) THEN type = TYP_DCMPLX
   IF (KEYWORD_SET(double) EQ true) THEN begin
      IF (type EQ TYP_FLOAT) THEN type = TYP_DOUBLE
      IF (type EQ TYP_COMPLEX) THEN type = TYP_DCMPLX
   END
   ;
   ; Setup the parameters for the call to the system function.
   ;
   ; Floating point arguments and keywords

   ; Get space for the vectors.
   ; The variable scv is also used later to determine how to post process
   ; the vectors.
   IF (ARG_PRESENT(vectors)) THEN BEGIN 
      scv = IMSL_2
      IF ((type EQ typ_float) OR (type EQ typ_complex)) $
           THEN vectors_spc = fltarr(scv*n, n) ELSE vectors_spc = dblarr(scv*n, n)
   END 
   ;
   ; Get the space for the results.  
   ;
   ; ALPHA is complex in all cases.
   sc_alpha = IMSL_2
   IF ((type EQ TYP_FLOAT) OR (type EQ  TYP_COMPLEX)) THEN BEGIN 
      alpha_spc = fltarr(sc_alpha*n)
   END ELSE BEGIN
      alpha_spc = dblarr(sc_alpha*n)
   END
   ; BETA is complex only if type si TYP_COMPLEX or TYP_DCMPLX.
   sc_beta = 1
   IF ((type EQ TYP_DCMPLX) OR (type EQ  TYP_COMPLEX)) THEN  sc_beta = IMSL_2
   IF ((type EQ TYP_FLOAT) OR (type EQ  TYP_COMPLEX)) THEN BEGIN 
      beta_spc = fltarr(sc_beta*n)
   END ELSE BEGIN
      beta_spc = dblarr(sc_beta*n)
   END
   ;
   ; Convert A and B
   ;
   a_cvt = imsl_cvt_arr(a, type)
   b_cvt = imsl_cvt_arr(b, type)
   ;
   ; Call the system function.
   err_status = 0L 
   MATHSTAT_222, type, err_status, a_cvt, b_cvt,  n, $ 
                              alpha_spc, $
                              beta_spc, $
                              vectors_spc
   ;
   ; Now copy over all output keywords results.
   ;
   IF (ARG_PRESENT(vectors)) THEN  BEGIN 
       vectors = imsl_cvt_arr(vectors_spc, /back)
    END

   alpha = imsl_cvt_arr(alpha_spc, /back)
   IF ((type EQ TYP_DCMPLX) OR (type EQ  TYP_COMPLEX)) THEN $
     beta    =    imsl_cvt_arr(beta_spc,    /back) ELSE $
     beta    =    beta_spc
     
   ;
   ; return
   RETURN
END

                   
                   
                   

  
      

  
