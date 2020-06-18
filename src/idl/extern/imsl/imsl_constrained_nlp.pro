; $Id: //depot/idl/IDL_70/idl_src/libs/imsl/imsl_6.0/lib/imsl_constrained_nlp.pro#1 $
;
; Copyright (c) 1970-2006, VISUAL NUMERICS Inc. All Rights Reserved.
; This software is confidential information which is proprietary to and a
; trade secret of Visual Numerics Inc.  Use, duplication or disclosure is
; subject to the terms of an appropriate license agreement.
;
FUNCTION imsl_constrained_nlp, f, $                ;INPUT Scalar STRING
                m, $                          ;INPUT Scalar LONG
                n, $                          ;INPUT Scalar LONG
                itmax=itmax, $                ;INPUT Scalar LONG 
                xguess=xguess, $              ;INPUT 1-D array floating point
                xscale=xscale, $              ;INPUT 1-D array floating point
                difftype=difftype, $          ;INPUT Scalar LONG
                tau0=tau0, $                  ;INPUT Scalar floating point
                del0=del0, $                  ;INPUT Scalar floating point
                smallw=smallw, $              ;INPUT Scalar floating point
                delmin=delmin, $              ;INPUT Scalar floating point
                scfmax=scfmax, $              ;INPUT Scalar floating point
                epsdif=epsdif, $              ;INPUT Scalar floating point
                epsfcn=epsfcn, $              ;INPUT Scalar floating point
                taubnd=taubnd, $              ;INPUT Scalar floating point
                grad=grad, $                  ;INPUT Scalar STRING
                ibtype=ibtype, $              ;INPUT Scalar LONG 
                meq=meq, $                    ;INPUT Scalar LONG 
                obj=obj, $                    ;OUTPUT Scalar floating point
                xlb=xlb, $                    ;IN/OUTPUT 1-D array floating point
                xub=xub, $                    ;IN/OUTPUT 1-D array floating point
                double = double               ;INPUT Scalar ON/OFF flag

@imsl_init.pro
   ON_ERROR, on_err_action
   ;
   ; Error checking.
   ; - F must be a scalar string.
   ; - if XGUESS is supplied, it must be a 1-D array of length n.
   ; - if XSCALE is supplied, it must be a 1-D array of length n.
   ; - if XLB is input (ibtype == 0 or 3) then it must be a 1-D array
   ;   of length n.
   ; - if XUB is input (ibtype == 0 or 3) then it must be a 1-D array
   ;   of length n.
   ; - if MEQ is not supplied, then set meq = m (Default)
   ;                       
   ;                       
   ; We do the following for both XLB and XUB:
   ; - if input:( ibtype eq 0 or 3)
   ;     o   if not present, then get a temporary array and fill it with
   ;         defaults.
   ;     o   if present, convert the array to the correct type, 
   ;         and use in the call to C/Math.
   ; - else  (Output only: ibtype eq 1 or 2)
   ;     o   get space to use in call to C/Math.  Since Ouput only,
   ;         the contents do no matter.
   ;
   ;                       
   nargs = n_params()
   IF (nargs NE 3) THEN message, 'Incorrect number of arguments.'
   size_f = IMSL_SIZE(f)
   IF ((N_ELEMENTS(f) NE 1) OR (size_f(N_ELEMENTS(size_f)-2) NE 7)) THEN $
     message, 'F must be a scalar string.'
   m_cvt = IMSL_LONG(m(0))
   n_cvt = IMSL_LONG(n(0))
   IF (n_cvt LT 1) THEN message, 'The number of functions must be positive.'
   n_cvt = IMSL_LONG(n(0))
   IF (n_cvt LT 1) THEN message, 'The number of variables must be positive.'
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
   IF ((KEYWORD_SET(meq) + ARG_PRESENT(meq)) GE 1) THEN meq_cvt = IMSL_LONG(meq(0)) ELSE meq_cvt = m_cvt
   IF (KEYWORD_SET(grad)) THEN BEGIN
      size_grad = IMSL_SIZE(grad)
      IF ((N_ELEMENTS(grad) NE 1) OR (size_grad(N_ELEMENTS(size_grad)-2) NE 7)) THEN $
        message, 'GRAD must be a scalar string.'
   END
   ; 
   ; Decide on what precision to use.
   type = TYP_FLOAT
   IF (KEYWORD_SET(double) EQ true) THEN type = TYP_DOUBLE
   ;
   ; Setup the parameters for the call to the system function.
   ; Input LONG keyword(s)
   ;
   IF (KEYWORD_SET(difftype)) THEN difftype_cvt = IMSL_LONG(difftype(0))
   IF (KEYWORD_SET(itmax)) THEN itmax_cvt = IMSL_LONG(itmax(0))
   IF (KEYWORD_SET(ibtype)) THEN ibtype_cvt = IMSL_LONG(ibtype(0)) ELSE ibtype_cvt = IMSL_0
   ;
   ; Floating point arguments and keywords   
   IF (type EQ TYP_DOUBLE) THEN BEGIN
      result = dblarr(n_cvt)
      IF (KEYWORD_SET(xguess)) THEN xguess_cvt = double(xguess)
      IF (KEYWORD_SET(xscale)) THEN xscale_cvt = double(xscale)
      IF (KEYWORD_SET(tau0)) THEN tau0_cvt = double(tau0(0))
      IF (KEYWORD_SET(del0)) THEN del0_cvt = double(del0(0))
      IF (KEYWORD_SET(smallw)) THEN smallw_cvt = double(smallw(0))
      IF (KEYWORD_SET(scfmax)) THEN scfmax_cvt = double(scfmax(0))
      IF (KEYWORD_SET(delmin)) THEN delmin_cvt = double(delmin(0))
      IF (KEYWORD_SET(epsdif)) THEN epsdif_cvt = double(epsdif(0))
      IF (KEYWORD_SET(epsfcn)) THEN epsfcn_cvt = double(epsfcn(0))
      IF (KEYWORD_SET(taubnd)) THEN taubnd_cvt = double(taubnd(0))
      IF (ARG_PRESENT(obj)) THEN obj_spc = double(0.0)
      IF ((ibtype_cvt EQ 0) OR (ibtype_cvt EQ 3)) THEN BEGIN
         IF (KEYWORD_SET(xlb)) THEN BEGIN
            size_xlb = IMSL_SIZE(xlb)
            IF (size_xlb(0) NE 1) THEN message, 'XLB must be a 1-D array.'
            IF (size_xlb(1) NE n_cvt) THEN message, 'XLB is not the correct size.'
            xlb_cvt = double(xlb)
         END ELSE BEGIN
            xlb_cvt = dblarr(n_cvt)
            xlb_cvt(*) = -1.0e6
         END
         IF (KEYWORD_SET(xub)) THEN BEGIN
            size_xub = IMSL_SIZE(xub)
            IF (size_xub(0) NE 1) THEN message, 'XUB must be a 1-D array.'
            IF (size_xub(1) NE n_cvt) THEN message, 'XUB is not the correct size.'
            xub_cvt = double(xub)
         END ELSE BEGIN
            xub_cvt = dblarr(n_cvt)
            xub_cvt(*) = 1.0e6
         END
      END ELSE BEGIN            ;Output only
         xlb_cvt = dblarr(n)
         xub_cvt = dblarr(n)
      END
   END ELSE BEGIN
      result = fltarr(n_cvt)
      IF (KEYWORD_SET(xguess)) THEN xguess_cvt = float(xguess)
      IF (KEYWORD_SET(xscale)) THEN xscale_cvt = float(xscale)
      IF (KEYWORD_SET(tau0)) THEN tau0_cvt = float(tau0(0))
      IF (KEYWORD_SET(del0)) THEN del0_cvt = float(del0(0))
      IF (KEYWORD_SET(smallw)) THEN smallw_cvt = float(smallw(0))
      IF (KEYWORD_SET(delmin)) THEN delmin_cvt = float(delmin(0))
      IF (KEYWORD_SET(scfmax)) THEN scfmax_cvt = float(scfmax(0))
      IF (KEYWORD_SET(epsdif)) THEN epsdif_cvt = float(epsdif(0))
      IF (KEYWORD_SET(epsfcn)) THEN epsfcn_cvt = float(epsfcn(0))
      IF (KEYWORD_SET(taubnd)) THEN taubnd_cvt = float(taubnd(0))
      IF (ARG_PRESENT(obj)) THEN obj_spc = float(0.0)
      IF ((ibtype_cvt EQ 0) OR (ibtype_cvt EQ 3)) THEN BEGIN
         IF (KEYWORD_SET(xlb)) THEN BEGIN
            size_xlb = IMSL_SIZE(xlb)
            IF (size_xlb(0) NE 1) THEN message, 'XLB must be a 1-D array.'
            IF (size_xlb(1) NE n_cvt) THEN message, 'XLB is not the correct size.'
            xlb_cvt = float(xlb)
         END ELSE BEGIN
            xlb_cvt = fltarr(n_cvt)
            xlb_cvt(*) = -1.0e6
         END
         IF (KEYWORD_SET(xub)) THEN BEGIN
            size_xub = IMSL_SIZE(xub)
            IF (size_xub(0) NE 1) THEN message, 'XUB must be a 1-D array.'
            IF (size_xub(1) NE n_cvt) THEN message, 'XUB is not the correct size.'
            xub_cvt = float(xub)
         END ELSE BEGIN
            xub_cvt = fltarr(n_cvt)
            xub_cvt(*) = 1.0e6
         END
      END ELSE BEGIN            ;Output only
         xlb_cvt = fltarr(n)
         xub_cvt = fltarr(n)
      END
   END
   ;
   ; Call the system function.
   ;
   err_status = 0L

   MATHSTAT_501, type, err_status, $
                         f, $
                         m_cvt, $
                         n_cvt, $
                         itmax_cvt, $
                         xguess_cvt, $
                         xscale_cvt, $
                         difftype_cvt, $
                         tau0_cvt, $
                         del0_cvt, $
                         smallw_cvt, $
                         delmin_cvt, $
                         scfmax_cvt, $
                         epsdif_cvt, $
                         epsfcn_cvt, $
                         taubnd_cvt, $
                         grad, $
                         ibtype_cvt, $
                         meq_cvt, $
                         obj_spc, $
                         xlb_cvt, $
                         xub_cvt, $
                         result

   IF (ibtype_cvt GT 0) THEN BEGIN
      xlb = xlb_cvt
      xub = xub_cvt
   END
   IF (ARG_PRESENT(obj)) THEN obj = obj_spc
 RETURN, result
END

