; $Id: //depot/idl/IDL_70/idl_src/libs/imsl/imsl_6.0/lib/imsl_sp_gmres.pro#1 $
;
; Copyright (c) 1970-2006, VISUAL NUMERICS Inc. All Rights Reserved.
; This software is confidential information which is proprietary to and a
; trade secret of Visual Numerics Inc.  Use, duplication or disclosure is
; subject to the terms of an appropriate license agreement.
;
FUNCTION imsl_sp_gmres, amultp, $                      ;INPUT Scalar STRING
                   b, $                           ;INPUT 1-D array: floating point
                   itmax=itmax, $                 ;INPUT/OUTPUT Scalar LONG 
                   tolerance=tolerance,  $        ;INPUT Scalar floating point
                   precond=precond, $             ;INPUT Scalar STRING
                   max_krylov=max_krylov, $       ;INPUT Scalar LONG
                   double=double, $               ;INPUT Scalar ON/OFF Flag
                   hh_reorth=hh_reorth            ;INPUT Scalar LONG

@imsl_init.pro
   ON_ERROR, on_err_action
   ;
   ; Error checking.
   ; Checks on positional args:
   ; - AMULTP be a scalar string.
   ; - b must be a 1-D array. n = N_ELEMENTS(b)
   ;  
   ; If PRECOND is supplied, it must be a scalar string.
   
   ; GENERAL CHECKS:
   ; The following checks are made regardless of the number of positional
   ; arguments.
   ;
   nargs = n_params()
   IF (nargs NE 2)  THEN $
     message, 'Incorrect number of arguments.'     
   size_amultp = IMSL_SIZE(amultp)
   IF ((N_ELEMENTS(amultp) NE 1) OR (size_amultp(N_ELEMENTS(size_amultp)-2) NE 7)) THEN $
     message, 'AMULTP must be a scalar string.'

   size_b = IMSL_SIZE(b)
   IF (size_b(0) NE 1) THEN  $
     message, 'B must be a 1-D array.'
   n = IMSL_LONG(size_b(3))

   IF (KEYWORD_SET(precond)) THEN BEGIN
      size_precond = IMSL_SIZE(precond)
      IF ((N_ELEMENTS(precond) NE 1) OR $
          (size_precond(N_ELEMENTS(size_precond)-2) NE 7)) THEN $
        message, 'PRECOND must be a scalar string.'
   END
   
   ;
   ; Decide on what precision to use.
   type = TYP_FLOAT
   IF (size_b(N_ELEMENTS(size_b)-2) EQ  TYP_DOUBLE) THEN type = TYP_DOUBLE
   IF (KEYWORD_SET(double) EQ true) THEN type = TYP_DOUBLE
   ;
   ; Setup the parameters for the call to the system function.
   ;
   ; Input LONG argument(s)
   IF (KEYWORD_SET(hh_reorth) EQ TRUE) THEN $
     hh_reorth_cvt = IMSL_1
   IF (KEYWORD_SET(itmax) EQ TRUE) THEN $
     itmax_cvt = IMSL_LONG(itmax(0))
   IF (KEYWORD_SET(max_krylov) EQ TRUE) THEN $ 
     max_krylov_cvt = IMSL_LONG(max_krylov)
   ;
   ; Floating point arguments and keywords
   ;
   ; Result
   ;
   IF (type EQ TYP_DOUBLE) THEN BEGIN
      result = dblarr(n)
      ; Input 
      b_cvt = double(b)
      IF (KEYWORD_SET(tolerance)) THEN tolerance_cvt = double(tolerance(0))
   END ELSE BEGIN
      result = fltarr(n)
      ; Input 
      b_cvt = float(b)
      IF (KEYWORD_SET(tolerance)) THEN tolerance_cvt = float(tolerance(0))
   END
   
   ;
   ; Call the system function.
   ;
   err_status = 0L
   MATHSTAT_206,  type, err_status, amultp, b_cvt, n, $
                   itmax_cvt, $
                   tolerance_cvt, $
                   precond_cvt, $
                   max_krylov_cvt, $
                   hh_reorth_cvt, $
                   result

   ; Now copy over all output keywords results.
   IF (KEYWORD_SET(itmax)) THEN itmax = itmax_cvt
   ;
   ; Return.
   ;
   RETURN, result
END

                   
                   
                   

  
      

  
