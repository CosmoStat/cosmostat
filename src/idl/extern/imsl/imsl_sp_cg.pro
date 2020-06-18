; $Id: //depot/idl/IDL_70/idl_src/libs/imsl/imsl_6.0/lib/imsl_sp_cg.pro#1 $
;
; Copyright (c) 1970-2006, VISUAL NUMERICS Inc. All Rights Reserved.
; This software is confidential information which is proprietary to and a
; trade secret of Visual Numerics Inc.  Use, duplication or disclosure is
; subject to the terms of an appropriate license agreement.
;
FUNCTION imsl_sp_cg, amultp, $                         ;INPUT Scalar STRING
                   b, $                           ;INPUT 1-D array: floating point
                   itmax=itmax, $                 ;INPUT/OUTPUT Scalar LONG 
                   precond=precond, $             ;INPUT Scalar STRING
                   jacobi=jacobi, $               ;INPUT 1-D array LONG
                   double=DOUBLE, $                ;INPUT Scalar ON/OFF Flag
                   rel_err=rel_err                ;INPUT/OUTPUT Scalar floating point

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
   ; If JACOBI is supplied, it must me a 1-D array of length n.
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
   
   IF (KEYWORD_SET(jacobi)) THEN BEGIN 
      IF ((size_b(0) NE 1) OR (N_ELEMENTS(jacobi) NE n)) THEN  $
        message, 'JACOBI is not the correct size.'
   END
   ;
   ; Decide on what precision to use.
   ; Always use double precision.
   ;
   ; Decide on what precision to use.
   type = TYP_FLOAT
   IF (size_b(N_ELEMENTS(size_b)-2) EQ  TYP_DOUBLE) THEN type = TYP_DOUBLE
   IF (KEYWORD_SET(double) EQ true) THEN type = TYP_DOUBLE
   ;
   ; Setup the parameters for the call to the system function.
   ;
   ; Input LONG argument(s)
   IF (KEYWORD_SET(itmax) EQ TRUE) THEN $
     itmax_cvt = IMSL_LONG(itmax(0))
   ;
   ; Double arguments and keywords
   ;
   ; Result
   ;
   IF (type EQ TYP_DOUBLE) THEN BEGIN
      result = dblarr(n)
      ; 
      ; Input 
      b_cvt = double(b)
      IF (KEYWORD_SET(rel_err)) THEN rel_err_cvt = double(rel_err)
      IF (KEYWORD_SET(jacobi)) THEN jacobi_cvt = double(jacobi)
      ; Output
      IF (ARG_PRESENT(rel_err)) THEN rel_err_spc = 0.0d0
   END ELSE BEGIN
      result = fltarr(n)
      ; 
      ; Input 
      b_cvt = float(b)
      IF (KEYWORD_SET(rel_err)) THEN rel_err_cvt = float(rel_err)
      IF (KEYWORD_SET(jacobi)) THEN jacobi_cvt = float(jacobi)
      ; Output
      IF (ARG_PRESENT(rel_err)) THEN rel_err_spc = FLOAT(0.0)
   END
   
   ;
   ; Call the system function.
   ;
   err_status = 0L
   MATHSTAT_205,  type, err_status, amultp, b_cvt, n, $
                   itmax_cvt, $
                   rel_err_cvt, $
                   precond_cvt, $
                   jacobi_cvt, $
                   result

   ; Now copy over all output keywords results.
   IF (KEYWORD_SET(itmax)) THEN itmax = itmax_cvt
   IF (KEYWORD_SET(rel_err)) THEN rel_err = rel_err_cvt
   ;
   ; Return.
   ;
   RETURN, result
END

                   
                   
                   

  
      

  
