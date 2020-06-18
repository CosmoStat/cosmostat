; $Id: //depot/idl/IDL_70/idl_src/libs/imsl/imsl_6.0/lib/imsl_nlinlsq.pro#1 $
;
; Copyright (c) 1970-2006, VISUAL NUMERICS Inc. All Rights Reserved.
; This software is confidential information which is proprietary to and a
; trade secret of Visual Numerics Inc.  Use, duplication or disclosure is
; subject to the terms of an appropriate license agreement.
;
FUNCTION imsl_nlinlsq, f, $                    ;INPUT Scalar STRING
                m, $                          ;INPUT Scalar LONG
                n, $                          ;INPUT Scalar LONG
                xlb,  $                       ;INPUT 1-D array floating point
                xub,  $                       ;INPUT 1-D array floating point
                jacobian=jacobian, $          ;INPUT Scalar STRING
                double = double, $            ;INPUT Scalar ON/OFF flag
                intern_scale=intern_scale, $  ;INPUT Scalar ON/OFF flag
                max_evals=max_evals, $        ;INPUT Scalar LONG 
                itmax=itmax, $                ;INPUT Scalar LONG 
                n_digits=n_digits, $          ;INPUT Scalar LONG 
                max_jacobian=max_jacobian, $  ;INPUT Scalar LONG 
                max_step=max_step, $          ;INPUT Scalar floating point
                tol_afcn=tol_afcn, $          ;INPUT Scalar floating point
                tol_grad=tol_grad, $          ;INPUT Scalar floating point
                tol_rfcn=tol_rfcn, $          ;INPUT Scalar floating point
                tol_step=tol_step, $          ;INPUT Scalar floating point
                xguess=xguess, $              ;INPUT Scalar floating point
                tolerance=tolerance, $        ;INPUT Scalar floating point
                trust_region=trust_region, $  ;INPUT Scalar floating point
                fscale=fscale, $              ;INPUT 1-D array floating point
                xscale=xscale, $              ;INPUT 1-D array floating point
                rank=rank, $                  ;OUTPUT Scalar LONG 
                fjac=fjac, $                  ;OUTPUT 2-D array floating point
                fvec=fvec, $                  ;OUTPUT 1-D array floating point
                jtj_inverse=jtj_inverse       ;OUTPUT 2-D array floating point
   
@imsl_init.pro
   ON_ERROR, on_err_action
   ;
   ; Error checking.
   ; - F must be a scalar string.
   ; - The following keyword(s) are only valid if JACOBIAN is also set.
   ;    o MAX_JACOBIAN
   ; - if XGUESS is supplied, it must be a 1-D array of length n.
   ; - if XSCALE is supplied, it must be a 1-D array of length n.
   ; - if FSCALE is supplied, it must be a 1-D array of length m.
   ; - if XLB and XUB are supplied: then jtj_inverse is not allowed.
   ;       - then jtj_inverse is not allowed,
   ;       - XLB and XUB must  be a 1-D arrays of length n.
   ; - 
   ;                       
   nargs = n_params()
   IF ((nargs NE 3) AND (nargs NE 5)) THEN message, 'Incorrect number of arguments.'
   IF (nargs EQ 5) THEN bounded = 1 ELSE bounded = 0
   size_f = IMSL_SIZE(f)
   IF ((N_ELEMENTS(f) NE 1) OR (size_f(N_ELEMENTS(size_f)-2) NE 7)) THEN $
     message, 'F must be a scalar string.'
   m_cvt = IMSL_LONG(m(0))
   n_cvt = IMSL_LONG(n(0))
   IF (KEYWORD_SET(jacobian)) THEN BEGIN
      size_jacobian = IMSL_SIZE(jacobian)
      IF ((N_ELEMENTS(jacobian) NE 1) OR (size_jacobian(N_ELEMENTS(size_jacobian)-2) NE 7)) THEN $
        message, 'JACOBIAN must be a scalar string.'
   END ELSE BEGIN
      IF (KEYWORD_SET(max_jacobian)) THEN message, $
        'MAX_JACOBIAN is only valid when JACOBIAN is also supplied.'
   END
   IF (KEYWORD_SET(xguess)) THEN BEGIN
      size_xguess = IMSL_SIZE(xguess)
      IF (size_xguess(0) NE 1) THEN message, 'XGUESS must be a 1-D array.'
      IF (size_xguess(1) NE n_cvt) THEN message, 'XGUESS is not the correct size.'
   END
   IF (KEYWORD_SET(xscale)) THEN BEGIN
      size_xscale = IMSL_SIZE(xscale)
      IF (size_xscale(0) NE 1) THEN message, 'XSCALE must be a 1-D array.'
      IF (size_xscale(1) NE n_cvt) THEN message, 'XSCALE is not the correct size.'
   END
   IF (KEYWORD_SET(fscale)) THEN BEGIN
      size_fscale = IMSL_SIZE(fscale)
      IF (size_fscale(0) NE 1) THEN message, 'FSCALE must be a 1-D array.'
      IF (size_fscale(1) NE m_cvt) THEN message, 'FSCALE is not the correct size.'
   END
   IF (bounded NE 0) THEN BEGIN 
      IF (ARG_PRESENT(jtj_inverse)) THEN $
          MESSAGE,   'JTJ_INVERSE is not allowed in a problem with bounds.'
      size_xlb = IMSL_SIZE(xlb)
      IF (size_xlb(0) NE 1) THEN message, 'XLB must be a 1-D array.'
      IF (N_ELEMENTS(xlb) NE n_cvt) THEN message, 'XLB is not the correct size'
      size_xub = IMSL_SIZE(xub)
      IF (size_xub(0) NE 1) THEN message, 'XUB must be a 1-D array.'
      IF (N_ELEMENTS(xub) NE n_cvt) THEN MESSAGE,  'XUB is not the correct size'
   END
   
   ; 
   ; Decide on what precision to use.
   type = TYP_FLOAT
   IF (KEYWORD_SET(double) EQ true) THEN type = TYP_DOUBLE
   IF (bounded NE 0) THEN BEGIN 
      IF (size_xlb(N_ELEMENTS(size_xlb)-2) EQ  TYP_DOUBLE) THEN type = TYP_DOUBLE
      IF (size_xub(N_ELEMENTS(size_xub)-2) EQ  TYP_DOUBLE) THEN type = TYP_DOUBLE
   END
   
   ;
   ; Setup the parameters for the call to the system function.
   ; Input LONG keyword(s)
   ;
   IF (KEYWORD_SET(itmax)) THEN itmax_cvt = IMSL_LONG(itmax(0))
   IF (KEYWORD_SET(max_evals)) THEN max_evals_cvt = IMSL_LONG(max_evals(0))
   IF (KEYWORD_SET(max_jacobian)) THEN max_jac_cvt = IMSL_LONG(max_jacobian(0))
   IF (KEYWORD_SET(n_digits)) THEN n_digits_cvt = IMSL_LONG(n_digits(0))
   IF (KEYWORD_SET(intern_scale)) THEN  intern_scale_cvt = IMSL_1
   IF (ARG_PRESENT(rank)) THEN  rank_spc = IMSL_0
   ;
   ; Floating point arguments and keywords   
   IF (type EQ TYP_DOUBLE) THEN BEGIN
      result = dblarr(n_cvt)
      IF (KEYWORD_SET(max_step)) THEN max_step_cvt = double(max_step(0))
      IF (KEYWORD_SET(tol_grad)) THEN tol_grad_cvt = double(tol_grad(0))
      IF (KEYWORD_SET(fscale)) THEN fscale_cvt = double(fscale)
      IF (KEYWORD_SET(tol_grad)) THEN tol_grad_cvt = double(tol_grad(0))
      IF (KEYWORD_SET(tol_afcn)) THEN tol_afcn_cvt = double(tol_afcn(0))
      IF (KEYWORD_SET(tol_rfcn)) THEN tol_rfcn_cvt = double(tol_rfcn(0))
      IF (KEYWORD_SET(tol_step)) THEN tol_step_cvt = double(tol_step(0))
      IF (KEYWORD_SET(tolerance)) THEN tolerance_cvt = double(tolerance(0))
      IF (KEYWORD_SET(trust_region)) THEN trust_cvt = double(trust_region(0))
      IF (KEYWORD_SET(xguess)) THEN xguess_cvt = double(xguess)
      IF (KEYWORD_SET(xscale)) THEN xscale_cvt = double(xscale)
      IF (ARG_PRESENT(fjac)) THEN fjac_spc = dblarr(n_cvt, m_cvt)
      IF (ARG_PRESENT(fvec)) THEN fvec_spc = dblarr(m_cvt)
      IF (ARG_PRESENT(jtj_inverse)) THEN jtj_inverse_spc = dblarr(n_cvt, n_cvt)
   END ELSE BEGIN
      result = fltarr(n_cvt)
      IF (KEYWORD_SET(max_step)) THEN max_step_cvt = float(max_step(0))
      IF (KEYWORD_SET(tol_grad)) THEN tol_grad_cvt = float(tol_grad(0))
      IF (KEYWORD_SET(fscale)) THEN fscale_cvt = float(fscale)
      IF (KEYWORD_SET(tol_grad)) THEN tol_grad_cvt = float(tol_grad(0))
      IF (KEYWORD_SET(tol_afcn)) THEN tol_afcn_cvt = float(tol_afcn(0))
      IF (KEYWORD_SET(tol_rfcn)) THEN tol_rfcn_cvt = float(tol_rfcn(0))
      IF (KEYWORD_SET(tol_step)) THEN tol_step_cvt = float(tol_step(0))
      IF (KEYWORD_SET(tolerance)) THEN tolerance_cvt = float(tolerance(0))
      IF (KEYWORD_SET(trust_region)) THEN trust_cvt = float(trust_region(0))
      IF (KEYWORD_SET(xguess)) THEN xguess_cvt = float(xguess)
      IF (KEYWORD_SET(xscale)) THEN xscale_cvt = float(xscale)
      IF (ARG_PRESENT(fjac)) THEN fjac_spc = fltarr(n_cvt, m_cvt)
      IF (ARG_PRESENT(fvec)) THEN fvec_spc = fltarr(m_cvt)
      IF (ARG_PRESENT(jtj_inverse)) THEN jtj_inverse_spc = fltarr(n_cvt, n_cvt)
   END
   IF (bounded NE 0) THEN BEGIN 
      xlb_cvt = imsl_cvt_arr(xlb, type)
      xub_cvt = imsl_cvt_arr(xub, type)
   END
   ;
   ; Call the system function.
   ;
   err_status = 0L
   MATHSTAT_163, type, err_status, $
                         f, $
                         m_cvt, $
                         n_cvt, $
                         jacobian, $
                         intern_scale_cvt, $
                         max_evals_cvt, $
                         itmax_cvt, $
                         n_digits_cvt, $
                         max_jac_cvt, $
                         max_step_cvt, $
                         tol_afcn_cvt, $
                         tol_grad_cvt, $
                         tol_rfcn_cvt, $
                         tol_step_cvt, $
                         xguess_cvt, $
                         tolerance_cvt, $
                         trust_cvt, $
                         fscale_cvt, $
                         xscale_cvt, $
                         rank_spc, $
                         fjac_spc, $
                         fvec_spc, $
                         jtj_inverse_spc, $
                         xlb_cvt, $
                         xub_cvt, $
                         result
   IF (ARG_PRESENT(rank)) THEN rank = rank_spc
   IF (ARG_PRESENT(fjac)) THEN fjac = fjac_spc
   IF (ARG_PRESENT(fvec)) THEN fvec = fvec_spc
   IF (ARG_PRESENT(jtj_inverse)) THEN jtj_inverse = jtj_inverse_spc
 RETURN, result
END

