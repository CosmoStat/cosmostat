; $Id: //depot/idl/IDL_70/idl_src/libs/imsl/imsl_6.0/lib/imsl_quadprog.pro#1 $      
;
; Copyright (c) 1970-2006, VISUAL NUMERICS Inc. All Rights Reserved.
; This software is confidential information which is proprietary to and a
; trade secret of Visual Numerics Inc.  Use, duplication or disclosure is
; subject to the terms of an appropriate license agreement.
;
FUNCTION imsl_quadprog, a, $                     ;INPUT 2-D array:floating point
                b, $                          ;INPUT 1-D array:floating point
                g, $                          ;INPUT 1-D array:floating point
                h, $                          ;INPUT 2-D array:floating point
                obj=obj, $                    ;OUTPUT Scalar floating point
                dual=dual, $                  ;OUTPUT 1-D array:floating point
                double=double, $              ;INPUT Scalar ON/OFF flag
                diag=diag, $                  ;OUTPUT Scalar floating point
                meq=meq                       ;INPUT Scalar LONG

@imsl_init.pro
   ON_ERROR, on_err_action
   ;
   ; Error checking:
   ;    - A must be 2-D.  (m x n)
   ;    - B must be 1-D of length m.
   ;    - G must be 1-D of length n.
   ;    - H must be 2-D.of size (n x n)
   ;    - If MEQ is not specified, then set the default, meq = m;
   ;          
   nargs = n_params()
   IF (nargs NE 4) THEN  message, 'Incorrect number of arguments.'
   size_a = IMSL_SIZE(a)
   IF (size_a(0) NE 2) THEN $
     message, 'A must be a 2-D array.'
   m = IMSL_LONG(size_a(1))
   n = IMSL_LONG(size_a(2))
   size_b = IMSL_SIZE(b)
   IF (size_b(0) NE 1) THEN $
     message, 'B must be a 1-D array.'
   IF (size_b(1) NE m) THEN $
     message, 'B is not the correct size.'
   size_g = IMSL_SIZE(g)
   IF (size_g(0) NE 1) THEN $
     message, 'G must be a 1-D array.'
   IF (size_g(1) NE n) THEN $
     message, 'G is not the correct size.'
   size_h = IMSL_SIZE(h)
   IF (size_h(0) NE 2) THEN $
     message, 'H must be a 2-D array.'
   IF ((size_h(1) NE n) OR (size_h(2) NE n)) THEN $
     message, 'H is not the correct size.'
   ;
   ; Decide on what precision to use.
   ;
   type = TYP_FLOAT
   IF (size_a(N_ELEMENTS(size_a)-2) EQ  TYP_DOUBLE) THEN type = TYP_DOUBLE
   IF (size_b(N_ELEMENTS(size_b)-2) EQ  TYP_DOUBLE) THEN type = TYP_DOUBLE
   IF (size_g(N_ELEMENTS(size_g)-2) EQ  TYP_DOUBLE) THEN type = TYP_DOUBLE
   IF (size_h(N_ELEMENTS(size_h)-2) EQ  TYP_DOUBLE) THEN type = TYP_DOUBLE
   IF (KEYWORD_SET(double) EQ true) THEN type = TYP_DOUBLE
   ;
   ; Setup the parameters for the call to the system function.
   ;
   if (arg_present(meq)) then meq_cvt = IMSL_LONG(meq(0)) ELSE meq_cvt = m
   ;
   ; Floating point arguments and keywords
   IF ((type eq TYP_DOUBLE) OR (type EQ TYP_DCMPLX)) THEN BEGIN
      result = dblarr(n)
      a_cvt = double(transpose(a))
      b_cvt = double(b)
      g_cvt = double(g)
      h_cvt = double(transpose(h))
      IF (ARG_PRESENT(obj)) THEN obj_spc = double(0.0)
      IF (ARG_PRESENT(dual)) THEN dual_spc = dblarr(n)
      IF (ARG_PRESENT(diag)) THEN diag_spc = double(0.0)
   END ELSE BEGIN
      result = fltarr(n)
      a_cvt = float(transpose(a))
      b_cvt = float(b)
      g_cvt = float(g)
      h_cvt = float(transpose(h))
      IF (ARG_PRESENT(obj)) THEN obj_spc = float(0.0)
      IF (ARG_PRESENT(dual)) THEN dual_spc = fltarr(n)
      IF (ARG_PRESENT(diag)) THEN diag_spc = float(0.0)
   END
   ;
   ; Call the system function.
   ;
   err_status = 0L
   MATHSTAT_178, type, err_status, a_cvt, $
                              b_cvt, $
                              g_cvt, $
                              h_cvt, $
                              m, $
                              n, $
                              obj_spc, $
                              dual_spc, $
                              diag_spc, $
                              meq_cvt, $
                              result
   
   IF (ARG_PRESENT(diag)) THEN diag = diag_spc
   IF (ARG_PRESENT(dual)) THEN dual = dual_spc
   IF (ARG_PRESENT(obj)) THEN obj = obj_spc
   ;
   ; Return
   RETURN, result
END
