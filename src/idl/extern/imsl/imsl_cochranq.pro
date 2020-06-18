; $Id: //depot/idl/IDL_70/idl_src/libs/imsl/imsl_6.0/lib/imsl_cochranq.pro#1 $
;
; Copyright (c) 1970-2006, VISUAL NUMERICS Inc. All Rights Reserved.
; This software is confidential information which is proprietary to and a
; trade secret of Visual Numerics Inc.  Use, duplication or disclosure is
; subject to the terms of an appropriate license agreement.
;
FUNCTION imsl_cochranq, x, $              ;INPUT 2-D array: floating point 
           q=q,  $                        ;OUTPUT Scalar floating point
           double=double                  ;INPUT Scalar ON/OFF flag
           
@imsl_init.pro
   ON_ERROR, on_err_action
   ;
   ; Error checking:
   ;  - X must be either a 1-D or 2-D array. (m x n)
   ;          
   nargs = n_params()
   IF (nargs NE 1)  THEN message, 'Incorrect number of arguments.'
     
   size_x = IMSL_LONG(size(x))
   IF ((size_x(0) NE 1) AND (size_x(0) NE 2)) THEN $
     message, 'X must be a 1-D or 2-D array.'
   m = IMSL_LONG(size_x(1))
   IF (size_x(0) EQ 2) THEN n = IMSL_LONG(size_x(2)) ELSE n = IMSL_1
   ;
   ; Decide on what precision to use.
   ;
   type = TYP_FLOAT
   IF (size_x(N_ELEMENTS(size_x)-2) EQ  TYP_DOUBLE) THEN type = TYP_DOUBLE
   IF (KEYWORD_SET(double) EQ true) THEN type = TYP_DOUBLE
   ;
   ; Setup the parameters for the call to the system function.
   ;
   ; Floating point arguments and keywords
   ;
   IF (type EQ TYP_DOUBLE) THEN BEGIN
      result = DOUBLE(0.)
      q_spc = double(0.)
   END ELSE BEGIN
      result = FLOAT(0.)
      q_spc = FLOAT(0.)
   END
   x_cvt = imsl_cvt_arr(x, type)
   ;
   ; Call the system function.
   ;
   err_status = 0L
   MATHSTAT_252,  type,  err_status,  x_cvt,  m,  n,  q_spc,  result
   ;
   q = q_spc
   ;
   ; Return
   RETURN, result
END
