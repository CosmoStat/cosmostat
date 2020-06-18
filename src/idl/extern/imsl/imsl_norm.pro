; $Id: //depot/idl/IDL_70/idl_src/libs/imsl/imsl_6.0/lib/imsl_norm.pro#1 $
;
; Copyright (c) 1970-2006, VISUAL NUMERICS Inc. All Rights Reserved.
; This software is confidential information which is proprietary to and a
; trade secret of Visual Numerics Inc.  Use, duplication or disclosure is
; subject to the terms of an appropriate license agreement.
;
FUNCTION imsl_norm, x, $                             ;INPUT 1-D array: floating point
               y, $                               ;INPUT 1-D array: floating point
               index_max=index_max, $             ;OUTPUT Scalar LONG
               inf=inf, $                         ;INPUT Scalar ON/OFF flag
               one=one, $                         ;INPUT Scalar ON/OFF flag
               double=double                      ;INPUT Scalar ON/OFF flag


@imsl_init.pro
   ON_ERROR, on_err_action

   ;
   ; Error checking.  
   ;    - x must be 1-D. Set n = n_elements(x)
   ;    - if y is present, then must be 1-D and of length (n). 
   ;    - ONE and INF are mutually exclusive.
   ;    - ONE and INDEX_MAX are mutually exclusive.
   ;    - If INDEX_MAX is used, then INF must also be present.
   ;
   nargs = n_params()
   IF ((nargs NE 1) and (nargs NE 2)) THEN $
         message, 'Incorrect number of arguments.'
   ;
   size_x = IMSL_SIZE(x)
   IF (size_x(0) NE 1) THEN BEGIN 
      message, 'X must be a 1-D array.'
   END
   n = size_x(1)
   IF (nargs EQ 2) THEN BEGIN
      size_y = IMSL_SIZE(y)
      IF ((N_ELEMENTS(y) NE n) OR (size_y(0) NE 1)) THEN BEGIN 
         message, 'Y must have the same dimensions as X.'
      END
   END
   IF ((KEYWORD_SET(one) + KEYWORD_SET(inf)) EQ 2) THEN $
     message, 'The keywords ONE and INF are mutually exclusive.'
   IF ((KEYWORD_SET(one) + ARG_PRESENT(index_max)) EQ 2) THEN $
     message, 'The keywords ONE and INDEX_MAX are mutually exclusive.'
   IF (ARG_PRESENT(index_max)) THEN $
     IF (NOT KEYWORD_SET(inf)) THEN $
     message, 'If INDEX_MAX is set, then INF must also be set.'
   ;
   ;ERROR CHECKING COMPLETE.
   ;
   ; Decide on what precision to use.
   ; Use the highest precision of the input argument(s).
   ;
   type = TYP_FLOAT
   IF (size_x(N_ELEMENTS(size_x)-2) GT TYP_FLOAT) THEN type = TYP_DOUBLE
   IF (nargs EQ 2) THEN $
     IF (size_y(N_ELEMENTS(size_y)-2) GT TYP_FLOAT) THEN type = TYP_DOUBLE
   IF (KEYWORD_SET(double) EQ true) THEN type = TYP_DOUBLE
   ;
   ; Setup the parameters for the call to the system function.
   ;
   ; Input LONG keyword(s)
   IF (KEYWORD_SET(one) EQ TRUE) THEN $
     one_cvt = IMSL_1
   IF (KEYWORD_SET(inf) EQ TRUE) THEN $
     inf_cvt = IMSL_1
   ;
   ; Input LONG keyword(s).
   ; Always ask for index_max, return to caller only when asked for.
   index_max_spc = IMSL_LONG(-1)
     
   ;
   ; Floating point arguments and keywords
   IF (type EQ TYP_DOUBLE) THEN BEGIN
      ; Result
      result = double(0.0)
      ; 
      ; Input
      x_cvt = double(x)
      IF (nargs EQ 2) THEN y_cvt = double(y)
   END ELSE BEGIN
      ; Result
      result = float(0.0)
      ; 
      ; Input
      x_cvt = float(x)
      IF (nargs EQ 2) THEN y_cvt = float(y)
   END
   ;
   ; Call the system function.
   ;
   err_status = 0L
   MATHSTAT_166,   type, err_status, $
               x_cvt, $
               y_cvt, $
               n, $
               one_cvt, $
               inf_cvt, $
               index_max_spc, $
               result
                 
   ;
   ; Now copy over all output keywords results.
   ;
      IF (ARG_PRESENT(index_max) EQ TRUE) THEN $
         index_max = index_max_spc
   ;
   ; Return.
   ;
   RETURN, result
END
   

