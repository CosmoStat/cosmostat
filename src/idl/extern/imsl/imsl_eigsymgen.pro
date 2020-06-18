; $Id: //depot/idl/IDL_70/idl_src/libs/imsl/imsl_6.0/lib/imsl_eigsymgen.pro#1 $   
;
; Copyright (c) 1970-2006, VISUAL NUMERICS Inc. All Rights Reserved.
; This software is confidential information which is proprietary to and a
; trade secret of Visual Numerics Inc.  Use, duplication or disclosure is
; subject to the terms of an appropriate license agreement.
;
   
FUNCTION imsl_eigsymgen, a, $                 ;INPUT 2-D array: floating point 
                b, $                       ;INPUT 2-D array: floating point 
                double=double, $           ;INPUT Scalar ON/OFF flag
                vectors=vectors            ;OUTPUT 2-D array: floating point 

@imsl_init.pro
   ON_ERROR, on_err_action
   ;
   ; Error checking:
   ;    - A must be 2-D.
   ;    - A must be square. Set n = p_argv[0]->value.arr->dim[0]
   ;    - B must be 2-D.
   ;    - B must be square and (n x n).
   ;
   nargs = n_params()
   IF (nargs NE 2) THEN message, 'Incorrect number of arguments.'
   size_a = IMSL_SIZE(a)
   IF (size_a(0) NE 2) THEN message, 'A must be a 2-D square array.'
   IF (size_a(1) NE size_a(2)) THEN message, 'A must be a 2-D square array.' $
     ELSE n = IMSL_LONG(size_a(1))
   size_b = IMSL_SIZE(b)
   IF (size_b(0) NE 2) THEN message, 'B must be a 2-D square array.'
   IF ((size_b(1) NE n) OR (size_b(2) NE n)) THEN $
     message, 'B must have the same dimensions as A.'
   ;
   ; Decide on what precision to use.
   ; Either TYP_FLOAT or TYP_DOUBLE *only*.
   type = TYP_FLOAT
   IF (size_a(N_ELEMENTS(size_a)-2) EQ  TYP_DOUBLE) THEN type = TYP_DOUBLE
   IF (size_b(N_ELEMENTS(size_b)-2) EQ  TYP_DOUBLE) THEN type = TYP_DOUBLE
   IF (KEYWORD_SET(double) EQ true) THEN type = TYP_DOUBLE
   ;
   ; Setup the parameters for the call to the system function.
   ;
   ; Input LONG keyword(s)
   ;
   ; Output LONG keyword(s)
   ;
   ; Floating point arguments and keywords
   ; Result vector.
   IF (type EQ TYP_DOUBLE) THEN BEGIN 
      result = dblarr(n)
      a_cvt = double(transpose(a))
      b_cvt = double(transpose(b))
      IF (ARG_PRESENT(vectors)) THEN vectors_spc = dblarr(n, n)
   END ELSE BEGIN
      result = fltarr(n)
      a_cvt = float(transpose(a))
      b_cvt = float(transpose(b))
      IF (ARG_PRESENT(vectors)) THEN vectors_spc = fltarr(n, n)
   END
   ;
   ; Convert A and B.
   ;
   a_cvt = imsl_cvt_arr(a, type)
   ;
   ; Call the system function.
   ;
   err_status = 0L
   MATHSTAT_134, type, err_status, a_cvt, b_cvt, n, $ 
                              vectors_spc, $
                              result
                           
   ;
   ; Now copy over all output keywords results.
   ;
   IF (ARG_PRESENT(vectors) EQ TRUE) THEN $
     vectors = transpose(vectors_spc)
   ;
   ; return
   RETURN, result
END

                   
                   
                   

  
      

  
