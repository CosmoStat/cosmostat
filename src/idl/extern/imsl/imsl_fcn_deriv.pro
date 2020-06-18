; $Id: //depot/idl/IDL_70/idl_src/libs/imsl/imsl_6.0/lib/imsl_fcn_deriv.pro#1 $
;
; Copyright (c) 1970-2006, VISUAL NUMERICS Inc. All Rights Reserved.
; This software is confidential information which is proprietary to and a
; trade secret of Visual Numerics Inc.  Use, duplication or disclosure is
; subject to the terms of an appropriate license agreement.
;
FUNCTION imsl_Fcn_deriv,   f,  $                 ;INPUT Scalar STRING
                      x,  $                 ;INPUT Scalar floating point
                      order=order,  $       ;INPUT Scalar LONG
                      stepsize=stepsize,  $ ;INPUT Scalar floating point
                      tolerance=tolerance, $ ;INPUT Scalar floating point
                      double=DOUBLE          ;INPUT Scalar ON/OFF flag

@imsl_init.pro
   ON_ERROR, on_err_action
   ;
   ; Error checking.
   ; - F must be a scalar string.
   ;                       
   nargs = n_params()
   IF (nargs NE 2) THEN message, 'Incorrect number of arguments.'
   size_f = IMSL_SIZE(f)
   IF ((N_ELEMENTS(f) NE 1) OR (size_f(N_ELEMENTS(size_f)-2) NE 7)) THEN $
     message, 'F must be a scalar string.'
   ; 
   ; Decide on what precision to use.
   type = TYP_FLOAT
   IF (KEYWORD_SET(DOUBLE) EQ true) THEN type  =  TYP_DOUBLE
   size_x = IMSL_SIZE(x)
   IF (size_x(N_ELEMENTS(size_x)-2) EQ  TYP_DOUBLE) THEN type = TYP_DOUBLE
   ;
   ; Setup the parameters for the call to the system function.
   ; Input LONG keyword(s)
   ;
   IF (KEYWORD_SET(order)) THEN order_cvt = IMSL_LONG(order(0)) ELSE order_cvt = IMSL_1
   ;
   ; Floating point arguments and keywords   
   IF (type EQ TYP_DOUBLE) THEN BEGIN
      mach = imsl_machine(/DOUBLE)
      result = DOUBLE(0.0)
      x_cvt = DOUBLE(x(0))
      IF (KEYWORD_SET(stepsize)) THEN stepsize_cvt  =  DOUBLE(stepsize(0)) $
        ELSE stepsize_cvt    =    DOUBLE(0.01)
      IF (KEYWORD_SET(tolerance)) THEN tolerance_cvt  =  DOUBLE(tolerance(0)) $
        ELSE tolerance_CVT = SQRT(SQRT(mach.max_rel_space))
   END ELSE BEGIN
      mach = imsl_machine(/FLOAT)
      result = 0.0
      x_cvt = FLOAT(x(0))
      IF (KEYWORD_SET(stepsize)) THEN stepsize_cvt = float(stepsize(0)) $
        ELSE stepsize_cvt    =    float(0.01)
      IF (KEYWORD_SET(tolerance)) THEN tolerance_cvt = float(tolerance(0)) $
        ELSE tolerance_cvt = SQRT(SQRT(mach.max_rel_space))
   END
   ;
   ; Call the system function.
   ;
   err_status = 0L  
   MATHSTAT_219, type, err_status, $
                              f, $
                              x_cvt, $
                              stepsize_cvt, $
                              tolerance_cvt, $
                              order_cvt, $
                              result
   
 RETURN, result
END

