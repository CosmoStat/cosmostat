; $Id: //depot/idl/IDL_70/idl_src/libs/imsl/imsl_6.0/lib/imsl_smoothdata1d.pro#1 $
;
; Copyright (c) 1970-2006, VISUAL NUMERICS Inc. All Rights Reserved.
; This software is confidential information which is proprietary to and a
; trade secret of Visual Numerics Inc.  Use, duplication or disclosure is
; subject to the terms of an appropriate license agreement.
;
FUNCTION imsl_smoothdata1d, x, $                       ;INPUT 1-D array: floating point
               y, $                               ;INPUT 1-D array: floating point
               itmax=itmax, $                     ;INPUT Scalar LONG
               distance=distance, $               ;INPUT Scalar floating point
               sc=sc, $                           ;INPUT Scalar floating point
               double=double                      ;INPUT Scalar ON/OFF flag

@imsl_init.pro
   ON_ERROR, on_err_action

   ;
   ; Error checking.  
   ;    - x must be 1-D. Set n = n_elements(x)
   ;    - y must be 1-D of length n.
   ;
   nargs = n_params()
   IF ((nargs NE 2)) THEN $
         message, 'Incorrect number of arguments.'
   ;
   size_x = IMSL_SIZE(x)
   IF (size_x(0) NE 1) THEN BEGIN 
      message, 'X must be a 1-D array.'
   END
   n = size_x(1)
   size_y = IMSL_SIZE(y)
   IF ((N_ELEMENTS(y) NE n) OR (size_y(0) NE 1)) THEN BEGIN 
         message, 'Y must have the same dimensions as X.'
   END
   ;
   ;ERROR CHECKING COMPLETE.
   ;
   ; Decide on what precision to use.
   ; Use the highest precision of the input argument(s).
   ;
   type = TYP_FLOAT
   IF (size_x(N_ELEMENTS(size_x)-2) GT TYP_FLOAT) THEN type = TYP_DOUBLE
   IF (size_y(N_ELEMENTS(size_y)-2) GT TYP_FLOAT) THEN type = TYP_DOUBLE
   IF (KEYWORD_SET(double) EQ true) THEN type = TYP_DOUBLE
   ;
   ; Setup the parameters for the call to the system function.
   ;
   ; Input LONG keyword(s)
   IF (KEYWORD_SET(itmax) EQ TRUE) THEN $
     itmax_cvt = IMSL_LONG(itmax(0))
   ;
   ; Floating point arguments and keywords
   IF (type EQ TYP_DOUBLE) THEN BEGIN
     IF (KEYWORD_SET(distance) EQ TRUE) THEN distance_cvt = double(distance(0))
     IF (KEYWORD_SET(sc) EQ TRUE) THEN sc_cvt = double(sc(0))
      ; Result
      result = dblarr(n)
      ; 
      ; Input
      x_cvt = double(x)
      y_cvt = double(y)
   END ELSE BEGIN
     IF (KEYWORD_SET(distance) EQ TRUE) THEN distance_cvt = float(distance(0))
     IF (KEYWORD_SET(sc) EQ TRUE) THEN sc_cvt = float(sc(0))
      ; Result
      result = fltarr(n)
      ; 
      ; Input
      x_cvt = float(x)
      y_cvt = float(y)
   END
   ;
   ; Call the system function.
   ;
   err_status = 0L
   MATHSTAT_214,   type, err_status, $
               x_cvt, $
               y_cvt, $
               n, $
               itmax_cvt, $
               distance_cvt, $
               sc_cvt, $
               result
                 
   ;
   ; Return.
   ;
   RETURN, result
END
   

