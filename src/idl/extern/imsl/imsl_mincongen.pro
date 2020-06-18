; $Id: //depot/idl/IDL_70/idl_src/libs/imsl/imsl_6.0/lib/imsl_mincongen.pro#1 $
;
; Copyright (c) 1970-2006, VISUAL NUMERICS Inc. All Rights Reserved.
; This software is confidential information which is proprietary to and a
; trade secret of Visual Numerics Inc.  Use, duplication or disclosure is
; subject to the terms of an appropriate license agreement.
;
FUNCTION imsl_Mincongen,   f,  $                   ;INPUT Scalar STRING
                      a,  $                   ;INPUT 2-D array floating point
                      b,  $                   ;INPUT 1-D array floating point
                      xlb,  $                 ;INPUT 1-D array floating point
                      xub,  $                 ;INPUT 1-D array floating point
                      meq=meq,  $             ;INPUT Scalar LONG 
                      xguess=xguess,  $       ;INPUT 1-D array floating point
                      max_fcn=max_fcn,  $     ;INPUT Scalar LONG 
                      grad=grad,  $           ;INPUT Scalar STRING
                      tolerance=tolerance,  $ ;INPUT Scalar floating point
                      obj=obj,   $            ;OUTPUT Scalar floating point
                      num_active=num_active,  $      ;OUTPUT Scalar LONG
                      active_const=active_const,  $  ;OUTPUT 1-D array LONG
                      lagrange_mult=lagrange_mult, $ ;OUTPUT 1-D array floating point
                      double  =  double              ;INPUT Scalar ON/OFF flag
   
@imsl_init.pro
   ON_ERROR, on_err_action
   ;
   ; Error checking.
   ; - F must be a scalar string.
   ; - A must be a 2-D array. set ncon and nvar based on A.
   ; - B must be a 1-D array of length ncon.
   ; - XLB must be a 1-D array of length nvar.
   ; - XUB must be a 1-D array of length nvar.
   ; 
   ; - if MEQ is not supplied, then set meq = 0 (Default)
   ; - if XGUESS is supplied, it must be a 1-D array of length nvar.
   ; - if GRAD is supplied, it must be a scalar string.
   ;                       
   nargs = n_params()
   IF (nargs NE 5) THEN message, 'Incorrect number of arguments.'
   size_f = IMSL_SIZE(f)
   IF ((N_ELEMENTS(f) NE 1) OR (size_f(N_ELEMENTS(size_f)-2) NE 7)) THEN $
     message, 'F must be a scalar string.'

   size_a = IMSL_SIZE(a)
   IF (size_a(0) NE 2) THEN message, 'A must be a 2-D array.'
   ncon = size_a(1)
   nvar = size_a(2)

   size_b = IMSL_SIZE(b)
   IF (size_b(0) NE 1) THEN message, 'B must be a 1-D array.'
   IF (N_ELEMENTS(b) NE ncon) THEN message, 'B is not the correct size'

   size_xlb = IMSL_SIZE(xlb)
   IF (size_xlb(0) NE 1) THEN message, 'XLB must be a 1-D array.'
   IF (N_ELEMENTS(xlb) NE nvar) THEN message, 'XLB is not the correct size'
   size_xub = IMSL_SIZE(xub)
   IF (size_xub(0) NE 1) THEN message, 'XUB must be a 1-D array.'
   IF (N_ELEMENTS(xub) NE nvar) THEN MESSAGE,  'XUB is not the correct size'

   IF (KEYWORD_SET(xguess)) THEN BEGIN
      size_xguess = IMSL_SIZE(xguess)
      IF (size_xguess(0) NE 1) THEN message, 'XGUESS must be a 1-D array.'
      IF (size_xguess(1) NE nvar) THEN message, 'XGUESS is not the correct size.'
   END

   IF (KEYWORD_SET(grad)) THEN BEGIN
      size_grad = IMSL_SIZE(grad)
      IF ((N_ELEMENTS(grad) NE 1) OR (size_grad(N_ELEMENTS(size_grad)-2) NE 7)) THEN $
        message, 'GRAD must be a scalar string.'
   END
   ; 
   ; Decide on what precision to use.
   type = TYP_FLOAT
   IF (size_a(N_ELEMENTS(size_a)-2) EQ  TYP_DOUBLE) THEN type = TYP_DOUBLE
   IF (size_b(N_ELEMENTS(size_b)-2) EQ  TYP_DOUBLE) THEN type = TYP_DOUBLE
   IF (size_xlb(N_ELEMENTS(size_xlb)-2) EQ  TYP_DOUBLE) THEN type = TYP_DOUBLE
   IF (size_xub(N_ELEMENTS(size_xub)-2) EQ  TYP_DOUBLE) THEN type = TYP_DOUBLE
   IF (KEYWORD_SET(double) EQ true) THEN type = TYP_DOUBLE
   ;
   ; Setup the parameters for the call to the system function.
   ; Input LONG keyword(s)
   ;
   num_active_spc = IMSL_0
   IF (KEYWORD_SET(meq)) THEN meq_cvt = IMSL_LONG(meq(0)) ELSE meq_cvt = IMSL_0
   IF (KEYWORD_SET(max_fcn)) THEN max_fcn_cvt = IMSL_LONG(max_fcn(0))
   IF (ARG_PRESENT(active_const)) THEN active_const_spc = IMSL_LONARR(ncon+2*nvar)
   ;
   ; Floating point arguments and keywords   
   IF (type EQ TYP_DOUBLE) THEN BEGIN
      result = dblarr(nvar)
      IF (KEYWORD_SET(xguess)) THEN xguess_cvt = double(xguess)
      IF (KEYWORD_SET(tolerance)) THEN tolerance_cvt = double(tolerance(0))
      IF (ARG_PRESENT(lagrange_mult)) THEN lagrange_m_spc = DBLARR(nvar)
      IF (ARG_PRESENT(obj)) THEN obj_spc = double(0.0)
   END ELSE BEGIN
      result = fltarr(nvar)
      IF (KEYWORD_SET(xguess)) THEN xguess_cvt = float(xguess)
      IF (KEYWORD_SET(tolerance)) THEN tolerance_cvt = float(tolerance(0))
      IF (ARG_PRESENT(lagrange_mult)) THEN lagrange_m_spc = FLTARR(nvar)
      IF (ARG_PRESENT(obj)) THEN obj_spc = float(0.0)
   END
   a_cvt = imsl_cvt_arr(a, type)
   b_cvt = imsl_cvt_arr(b, type)
   xlb_cvt = imsl_cvt_arr(xlb, type)
   xub_cvt = imsl_cvt_arr(xub, type)
   ;
   ; Call the system function.
   ;
   err_status = 0L
   MATHSTAT_218, type, err_status, $
                         f, $
                         nvar, $
                         ncon, $
                         a_cvt, $
                         b_cvt, $
                         xlb_cvt, $
                         xub_cvt, $
                         xguess_cvt, $
                         tolerance_cvt, $
                         grad, $
                         meq_cvt, $
                         obj_spc, $
                         num_active_spc, $
                         active_const_spc, $
                         lagrange_m_spc,  $
                         max_fcn_cvt, $
                         result
   num_active    =    num_active_spc
   IF (ARG_PRESENT(obj)) THEN obj  =  obj_spc
   IF (num_active GT 0) THEN BEGIN 
     IF (ARG_PRESENT(active_const)) THEN $
       active_const = active_const_spc(0:num_active-1)
     IF (ARG_PRESENT(lagrange_mult)) THEN $
       lagrange_mult = lagrange_m_spc(0:num_active-1)
   END
 RETURN, result
END

