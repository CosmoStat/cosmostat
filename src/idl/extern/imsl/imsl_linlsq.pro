; $Id: //depot/idl/IDL_70/idl_src/libs/imsl/imsl_6.0/lib/imsl_linlsq.pro#1 $
;
; Copyright (c) 1970-2006, VISUAL NUMERICS Inc. All Rights Reserved.
; This software is confidential information which is proprietary to and a
; trade secret of Visual Numerics Inc.  Use, duplication or disclosure is
; subject to the terms of an appropriate license agreement.
;
FUNCTION imsl_linlsq, b, $                       ;INPUT 1-D array: floating point 
                 a, $                       ;INPUT 2-D array: floating point 
                 c, $                       ;INPUT 2-D array: floating point 
                 bl, $                      ;INPUT 1-D array: floating point 
                 bu, $                      ;INPUT 1-D array: floating point 
                 contype, $                 ;INPUT 1-D array: LONG 
                 xlb=xlb, $                 ;INPUT 1-D array: floating point 
                 xub=xub, $                 ;INPUT 1-D array: floating point 
                 abs_tolerance = abs_tolerance, $ ;INPUT Scalar floating point
                 rel_tolerance=rel_tolerance, $   ;INPUT Scalar floating point
                 itmax=itmax , $            ;INPUT Scalar LONG
                 double=double, $           ;INPUT Scalar ON/OFF flag
                 residual=residual          ;OUTPUT 1-D array: floating point
@imsl_init.pro
   ON_ERROR, on_err_action
   ;
   ; Error checking:
   ;   - A must be a 2-D array of order (nra x nca)
   ;     nra and nca are set by this size of this argument.
   ;   - B must be a 1-D array of length nra
   ;   - C must be a 2-D array of order (ncon x nca)
   ;     ncon is set by this size of this argument.
   ;   - BL must be a 1-D array of length ncon.
   ;   - BU must be a 1-D array of length ncon.
   ;   - CONTYPE must be a 1-D array of length ncon.
   ;
   ;   If XLB/XUB are supplied:
   ;   - XLB/XUB must be a 1-D array of length nca.
   ;          
   nargs = n_params()
   IF (nargs NE 6)  THEN $
     message, 'Incorrect number of arguments.'
   size_a = IMSL_SIZE(a)
   IF (size_a(0) NE 2) THEN MESSAGE,  'A must be a 2-D array.'
   nra  =  IMSL_LONG(size_a(1))
   nca  =  IMSL_LONG(size_a(2))
   size_b  =  IMSL_LONG(SIZE(b))
   IF (size_b(0) NE 1) THEN MESSAGE,  'B must be a 1-D array.'
   IF (size_b(1) NE nra) THEN MESSAGE,  'B is not the correct length.'
   size_c = IMSL_SIZE(c)
   IF (size_c(0) NE 2) THEN MESSAGE,  'C must be a 2-D array.'
   ncon  =  IMSL_LONG(size_c(1))
   IF (size_c(2) NE nca) THEN MESSAGE, 'C is not the correct size.'
   size_bl  =  IMSL_LONG(SIZE(bl))
   IF (size_bl(0) NE 1) THEN MESSAGE,  'BL must be a 1-D array.'
   IF (size_bl(1) NE ncon) THEN MESSAGE,  'BL is not the correct length.'
   size_bu  =  IMSL_LONG(SIZE(bu))
   IF (size_bu(0) NE 1) THEN MESSAGE,  'BU must be a 1-D array.'
   IF (size_bu(1) NE ncon) THEN MESSAGE,  'BU is not the correct length.'
   size_contype  =  IMSL_LONG(SIZE(contype))
   IF (size_contype(0) NE 1) THEN MESSAGE,  'CONTYPE must be a 1-D array.'
   IF (size_contype(1) NE ncon) THEN MESSAGE,  'CONTYPE is not the correct length.'
   IF (KEYWORD_SET(xlb)) THEN BEGIN
      size_xlb   =   IMSL_LONG(SIZE(xlb))
      IF (size_xlb(0) NE 1) THEN MESSAGE,   'XLB must be a 1-D array.'
      IF (size_xlb(1) NE nca) THEN MESSAGE,    'XLB is not the correct length.'
   END ELSE BEGIN
      xlb     =     DBLARR(nca)
      xlb(*)  =  1.0d30
   END
   IF (KEYWORD_SET(xub)) THEN BEGIN
      size_xub   =   IMSL_LONG(SIZE(xub))
      IF (size_xub(0) NE 1) THEN MESSAGE,   'XUB must be a 1-D array.'
      IF (size_xub(1) NE nca) THEN MESSAGE,    'XUB is not the correct length.'
   END ELSE BEGIN
      xub     =     DBLARR(nca)
      xub(*)  =  -1.0d30
   END
   ;
   ; Decide on what precision to use.
   ;
   type = TYP_FLOAT
   IF (size_b(N_ELEMENTS(size_b)-2) EQ  TYP_DOUBLE) THEN type = TYP_DOUBLE
   IF (size_a(N_ELEMENTS(size_a)-2) EQ  TYP_DOUBLE) THEN type  =  TYP_DOUBLE
   IF (size_c(N_ELEMENTS(size_a)-2) EQ  TYP_DOUBLE) THEN type  =  TYP_DOUBLE
   IF (size_bl(N_ELEMENTS(size_a)-2) EQ  TYP_DOUBLE) THEN type  =  TYP_DOUBLE
   IF (size_bu(N_ELEMENTS(size_a)-2) EQ  TYP_DOUBLE) THEN type  =  TYP_DOUBLE
   IF (KEYWORD_SET(double) EQ true) THEN type = TYP_DOUBLE
   ;
   ; Setup the parameters for the call to the system function.
   ;
   ; Input LONG
   IF (KEYWORD_SET(itmax)) THEN itmax_cvt = IMSL_LONG(itmax) ELSE $
     itmax_cvt = 5L*IMSL_LONG(nra > nca)
   contype_cvt = IMSL_LONG(contype)
   ; Floating point arguments and keywords
   IF (type EQ TYP_DOUBLE) THEN BEGIN
      tmp = imsl_machine(/DOUBLE)
      result = dblarr(nca)
      abs_tol_cvt  =  SQRT(tmp.(3))
      rel_tol_cvt  =  SQRT(tmp.(3))
      IF (KEYWORD_SET(abs_tolerance)) THEN abs_tol_cvt = double(abs_tolerance(0))
      IF (KEYWORD_SET(rel_tolerance)) THEN rel_tol_cvt = double(rel_tolerance(0))
      residual_spc = dblarr(nra) ; Always get the residual.
    END ELSE BEGIN
      tmp = imsl_machine(/FLOAT)
      result = fltarr(nca)
      abs_tol_cvt  =  SQRT(tmp.(3))
      rel_tol_cvt  =  SQRT(tmp.(3))
      IF (KEYWORD_SET(abs_tolerance)) THEN abs_tol_cvt = float(abs_tolerance(0))
      IF (KEYWORD_SET(rel_tolerance)) THEN rel_tol_cvt = float(rel_tolerance(0))
      residual_spc = fltarr(nra) ; Always get the residual.
   END
   b_cvt   =   IMSL_CVT_ARR(b,  type)
   a_cvt   =   IMSL_CVT_ARR(a,  type)
   c_cvt   =   IMSL_CVT_ARR(c,  type)
   bl_cvt  =   IMSL_CVT_ARR(bl, type)
   bu_cvt  =   IMSL_CVT_ARR(bu, type)
   xlb_cvt =   IMSL_CVT_ARR(xlb,type)
   xub_cvt =   IMSL_CVT_ARR(xub,type)
   ;
   ; Call the system function.
   ;
   err_status = 0L
   MATHSTAT_221,   type, err_status,   $
                   nra, $
                   nca, $
                   ncon, $
                   b_cvt, $
                   a_cvt, $
                   c_cvt, $
                   bl_cvt, $
                   bu_cvt, $
                   xlb_cvt, $
                   xub_cvt, $
                   abs_tol_cvt, $
                   rel_tol_cvt, $
                   contype_cvt, $
                   residual_spc, $
                   itmax_cvt, $
                   result
   ;
   ; Now copy over all output keywords results.
   ;
   IF (ARG_PRESENT(residual)) THEN residual = residual_spc
   ;
   ; Return
   RETURN, result
END
