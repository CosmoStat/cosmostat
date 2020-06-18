; $Id: //depot/idl/IDL_70/idl_src/libs/imsl/imsl_6.0/lib/imsl_linprog.pro#1 $      
;
; Copyright (c) 1970-2006, VISUAL NUMERICS Inc. All Rights Reserved.
; This software is confidential information which is proprietary to and a
; trade secret of Visual Numerics Inc.  Use, duplication or disclosure is
; subject to the terms of an appropriate license agreement.
;
FUNCTION imsl_linprog, a, $                      ;INPUT 2-d array:floating point
                b, $                          ;INPUT 1-d array:floating point
                c, $                          ;INPUT 1-d array:floating point
                bu=bu, $                      ;INPUT 1-d array:floating point
                double=double, $              ;INPUT Scalar ON/OFF flag
                dual=dual, $                  ;OUTPUT 1-d array:floating point
                irtype=irtype, $              ;INPUT 1-d array:LONG
                itmax=itmax, $                ;INPUT Scalar LONG
                obj=obj, $                    ;OUTPUT Scalar floating point
                xlb=xlb, $                    ;INPUT 1-d array:floating point
                xub=xub                       ;INPUT 1-d array:floating point
 
@imsl_init.pro
   ON_ERROR, on_err_action
   ;
   ; Error checking:
   ;    - A must be 2-D.  (m x n)
   ;    - B must be 1-D of length m.
   ;    - C must be 1-D of length n.
   ;    - If BU is specified, it must be a 1-D array of length m.
   ;    - If IRTYPE is specified, it must be a 1-D array of length m.
   ;    - If XLB is specified, it must be a 1-D array of length n.
   ;    - If XUB is specified, it must be a 1-D array of length n.
   ;    - If ITMAX is not specified, then set the default, itmax = 10000;
   ;          
   nargs = n_params()
   IF (nargs NE 3) THEN  message, 'Incorrect number of arguments.'
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
   size_c = IMSL_SIZE(c)
   IF (size_c(0) NE 1) THEN $
     message, 'C must be a 1-D array.'
   IF (size_c(1) NE n) THEN $
     message, 'C is not the correct size.'
   IF (KEYWORD_SET(bu)) THEN BEGIN
      size_bu = IMSL_SIZE(bu)
      IF (size_bu(0) NE 1) THEN $
        message, 'BU must be a 1-D array.'
      IF (size_BU(1) NE m) THEN $
        message, 'BU is not the correct size.'
   END
   IF (KEYWORD_SET(irtype)) THEN BEGIN
      size_irtype = IMSL_SIZE(irtype)
      IF (size_irtype(0) NE 1) THEN $
        message, 'IRTYPE must be a 1-D array.'
      IF (size_IRTYPE(1) NE m) THEN $
        message, 'IRTYPE is not the correct size.'
   END
   IF (KEYWORD_SET(xlb)) THEN BEGIN
      size_xlb = IMSL_SIZE(xlb)
      IF (size_xlb(0) NE 1) THEN $
        message, 'XLB must be a 1-D array.'
      IF (size_XLB(1) NE n) THEN $
        message, 'XLB is not the correct size.'
   END
   IF (KEYWORD_SET(xub)) THEN BEGIN
      size_xub = IMSL_SIZE(xub)
      IF (size_xub(0) NE 1) THEN $
        message, 'XUB must be a 1-D array.'
      IF (size_XUB(1) NE n) THEN $
        message, 'XUB is not the correct size.'
   END
   ;
   ; Decide on what precision to use.
   ;
   type = TYP_FLOAT
   IF (size_a(N_ELEMENTS(size_a)-2) EQ  TYP_DOUBLE) THEN type = TYP_DOUBLE
   IF (size_b(N_ELEMENTS(size_b)-2) EQ  TYP_DOUBLE) THEN type = TYP_DOUBLE
   IF (size_c(N_ELEMENTS(size_c)-2) EQ  TYP_DOUBLE) THEN type = TYP_DOUBLE
   IF (KEYWORD_SET(double) EQ true) THEN type = TYP_DOUBLE
   ;
   ;
   ; Setup the parameters for the call to the system function.
   ;
   if (keyword_set(itmax)) then itmax_cvt = IMSL_LONG(itmax) ELSE itmax = IMSL_LONG(10000)
   if (keyword_set(irtype)) then irtype_cvt = IMSL_LONG(irtype) 
   ;
   ; Floating point arguments and keywords
   IF ((type eq TYP_DOUBLE) OR (type EQ TYP_DCMPLX)) THEN BEGIN
      result = dblarr(n)
      a_cvt = double(transpose(a))
      b_cvt = double(b)
      c_cvt = double(c)
      IF (KEYWORD_SET(bu)) THEN bu_cvt = double(bu)
      IF (KEYWORD_SET(xlb)) THEN xlb_cvt = double(xlb)
      IF (KEYWORD_SET(xub)) THEN xub_cvt = double(xub)
      IF (ARG_PRESENT(dual)) THEN dual_spc = dblarr(n)
      IF (ARG_PRESENT(obj)) THEN obj_spc = double(0.0)
   END ELSE BEGIN
      result = fltarr(n)
      a_cvt = float(transpose(a))
      b_cvt = float(b)
      c_cvt = float(c)
      IF (KEYWORD_SET(bu)) THEN bu_cvt = float(bu)
      IF (KEYWORD_SET(xlb)) THEN xlb_cvt = float(xlb)
      IF (KEYWORD_SET(xub)) THEN xub_cvt = float(xub)
      IF (ARG_PRESENT(dual)) THEN dual_spc = fltarr(n)
      IF (ARG_PRESENT(obj)) THEN obj_spc = float(0.0)
   END
   ;
   ; Call the system function.
   ;
   err_status = 0L
   MATHSTAT_154, type, err_status, a_cvt, $
                              b_cvt, $
                              c_cvt, $
                              m, $
                              n, $
                              bu_cvt, $
                              dual_spc, $
                              irtype_cvt, $
                              itmax_cvt, $
                              obj_spc, $
                              xlb_cvt, $
                              xub_cvt, $ 
                              result
   IF (ARG_PRESENT(dual)) THEN dual = dual_spc
   IF (ARG_PRESENT(obj)) THEN obj = obj_spc
   ;
   ; Return
   RETURN, result
END
