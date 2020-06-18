; $Id: //depot/idl/IDL_70/idl_src/libs/imsl/imsl_6.0/lib/imsl_eig.pro#1 $   
;
; Copyright (c) 1970-2006, VISUAL NUMERICS Inc. All Rights Reserved.
; This software is confidential information which is proprietary to and a
; trade secret of Visual Numerics Inc.  Use, duplication or disclosure is
; subject to the terms of an appropriate license agreement.
;
   
FUNCTION imsl_eig, a, $                       ;INPUT 2-D array: floating point 
                double=double, $           ;INPUT Scalar ON/OFF flag
                symmetric=symmetric, $     ;INPUT Scalar ON/OFF flag
                lower_limit=lower_limit, $ ;INPUT scalar floating point
                upper_limit=upper_limit, $ ;INPUT scalar floating point
                number=number, $           ;OUTPUT Scalar LONG
                vectors=vectors            ;OUTPUT 2-D array: floating point 
@imsl_init.pro
   ON_ERROR, on_err_action
   ;
   ; Error checking:
   ;    - A must be 2-D.
   ;    - A must be square
   ;    - The keywords LOWER_LIMIT and UPPER_LIMIT must be supplied
   ;      together, or not at all.
   ;    - The keywords LOWER_LIMIT and UPPER_LIMIT are only valid
   ;      if SYMMETRIC is also supplied.
   ;    - The keyword NUMBER is only valid if SYMMETRIC is 
   ;      also supplied.
   ;
   nargs = n_params()
   IF (nargs NE 1) THEN message, 'Incorrect number of arguments.'
   size_a = IMSL_SIZE(a)
   IF (size_a(0) NE 2) THEN message, 'A must be a 2-D square array.'
   IF (size_a(1) NE size_a(2)) THEN message, 'A must be a 2-D square array.' $
     ELSE n = IMSL_LONG(size_a(1))
   If ((KEYWORD_SET(lower_limit) + KEYWORD_SET(upper_limit)) EQ 1) THEN $
     message, 'The keywords LOWER_LIMIT and UPPER_LIMIT must be used together.'
   IF (KEYWORD_SET(symmetric) EQ FALSE) THEN BEGIN
      IF (KEYWORD_SET(lower_limit)) THEN $
        message, 'The keywords LOWER_LIMIT and UPPER_LIMIT are only valid ' + $
        'if the keyword SYMMETRIC is also supplied.'
      IF (KEYWORD_SET(number)) THEN $
        message, 'The keyword NUMBER only valid ' + $
        'if the keyword SYMMETRIC is also supplied.'
   END
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
   ; Input LONG keyword(s)
   IF (KEYWORD_SET(symmetric)) THEN symmetric_cvt = IMSL_1
   ;
   ; Output LONG keyword(s)
   IF (ARG_PRESENT(number)) THEN number_spc = -IMSL_1
   ;
   ; Floating point arguments and keywords

   ; Get space for the vectors.  If SYMMETRIC is set, then the vectors will be the
   ; the same precision as type.  Otherwise they will be either TYP_COMPLEX or
   ; TYP_DCMPLX.  The variable scv is also used later to determine how to post process
   ; the vectors.
   scv = IMSL_1
   IF (ARG_PRESENT(vectors)) THEN BEGIN 
      IF (KEYWORD_SET(symmetric))  THEN BEGIN
         IF ((type EQ TYP_COMPLEX) OR  (type EQ TYP_DCMPLX)) THEN scv = IMSL_2 ELSE scv = IMSL_1
      END ELSE scv = IMSL_2
      IF ((type EQ typ_float) OR (type EQ typ_complex)) $
           THEN vectors_spc = fltarr(scv*n, n) ELSE vectors_spc = dblarr(scv*n, n)
   END 
   ;
   ; Get the space for the result.  If SYMMETRIC is set, then the result is either
   ; TYP_FLOAT or TYP_DOUBLE.  Otherwise, the result is either TYP_COMPLEX or
   ; TYP_DCMPLX.
   ; We also convert the keywords LOWER_LIMIT and UPPER_LIMIT if they are present.
   ;
   IF (KEYWORD_SET(symmetric)) THEN sc = IMSL_1 ELSE sc = IMSL_2
   IF ((type EQ TYP_FLOAT) OR (type EQ  TYP_COMPLEX)) THEN BEGIN 
      result = fltarr(sc*n)
      IF (KEYWORD_SET(lower_limit)) THEN BEGIN
         lower_limit_cvt = float(lower_limit(0))
         upper_limit_cvt = float(upper_limit(0))
      END     
   END ELSE BEGIN
      result = dblarr(sc*n)
      IF (KEYWORD_SET(lower_limit)) THEN BEGIN
         lower_limit_cvt = double(lower_limit(0))
         upper_limit_cvt = double(upper_limit(0))
      END     
   END
   ;
   ; Convert A.
   ;
   a_cvt = imsl_cvt_arr(a, type)

   ;
   ; Call the system function.
   ;
   err_status = 0L
   MATHSTAT_133, type, err_status, a_cvt, n, $ 
                              symmetric_cvt, $
                              lower_limit_cvt, $
                              upper_limit_cvt, $
                              number_spc, $
                              vectors_spc, $
                              result
   ;
   ; Now copy over all output keywords results.
   ;
   IF (ARG_PRESENT(vectors)) THEN  BEGIN 
      IF (scv EQ 2) THEN vectors = imsl_cvt_arr(vectors_spc, /back) ELSE vectors = transpose(vectors_spc)
   END
   IF (ARG_PRESENT(number) EQ TRUE) THEN $
     number = number_spc
   ;
   ; return
   IF ((sc eq 2) and (not keyword_set(symmetric))) THEN RETURN, imsl_cvt_arr(result, /back) $
     ELSE RETURN, result
END

                   
                   
                   

  
      

  
