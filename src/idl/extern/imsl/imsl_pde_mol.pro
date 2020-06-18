; $Id: //depot/idl/IDL_70/idl_src/libs/imsl/imsl_6.0/lib/imsl_pde_mol.pro#1 $
;
; Copyright (c) 1970-2006, VISUAL NUMERICS Inc. All Rights Reserved.
; This software is confidential information which is proprietary to and a
; trade secret of Visual Numerics Inc.  Use, duplication or disclosure is
; subject to the terms of an appropriate license agreement.
;
FUNCTION imsl_pde_mol, t, $                     ;INPUT 1-D array: floating point 
                y, $                       ;INPUT 2-D array: floating point 
                xbreak, $                  ;INPUT 1-D array: floating point 
                f_ut, $                    ;INPUT Scalar STRING
                f_bc, $                    ;INPUT Scalar STRING
                tolerance=tolerance, $     ;INPUT Scalar floating point
                hinit=hinit, $             ;INPUT Scalar floating point
                deriv_init=deriv_init, $   ;INPUT 2-D array: floating point 
                double=double              ;INPUT Scalar ON/OFF flag

@imsl_init.pro
   ON_ERROR, on_err_action
   ;
   ; Error checking.
   ; - Y is a 2D array, (NPDE x NX)
   ; - t is a 1-D array. (n_t_vals = length.)
   ; - xbreak is a 1-D array of length nx.
   ; - f_ut must be a scalar string.
   ; - f_bc must be a scalar string.
   ; - If deriv_init is supplied, it must be a 2D array, (NPDE x NX)
   ;                       
   nargs = n_params()
   IF (nargs NE 5) THEN message, 'Incorrect number of arguments.'
   size_y = IMSL_SIZE(y)
   IF (size_y(0) NE 2) THEN message, 'Y must be a 2-D array.'
   npde = IMSL_LONG(size_y(1))
   nx = IMSL_LONG(size_y(2))
   size_t = IMSL_SIZE(t)
   size_xbreak = IMSL_SIZE(xbreak)
   IF (size_t(0) NE 1) THEN message, 'T must be a 1-D array.'
   n_t_vals = IMSL_LONG(size_t(1))
   IF (size_xbreak(0) NE 1) THEN message, 'XBREAK must be a 1-D array.'
  IF (N_ELEMENTS(xbreak) NE nx) THEN MESSAGE, 'XBREAK is not the correct size.'
   size_f_ut = IMSL_SIZE(f_ut)
   size_f_bc = IMSL_SIZE(f_bc)
   IF ((N_ELEMENTS(f_ut) NE 1) OR (size_f_ut(N_ELEMENTS(size_f_ut)-2) NE 7)) THEN $
     message, 'F_UT must be a scalar string.'
   IF ((N_ELEMENTS(f_bc) NE 1) OR (size_f_bc(N_ELEMENTS(size_f_bc)-2) NE 7)) THEN $
     message, 'F_BC must be a scalar string.'
   IF (KEYWORD_SET(deriv_init)) THEN BEGIN
      size_deriv_init = IMSL_SIZE(deriv_init)
      IF (size_deriv_init(0) NE 2) THEN message, 'DERIV_INIT must be a 2-D array.'
      IF ((size_deriv_init(1) NE npde) OR (size_deriv_init(2) NE nx)) THEN $
        MESSAGE, 'DERIVE_INIT is not the correct size.'
   END
   
   ; 
   ; Decide on what precision to use.
   type = TYP_FLOAT
   IF (size_t(N_ELEMENTS(size_t)-2) EQ  TYP_DOUBLE) THEN type = TYP_DOUBLE
   IF (size_xbreak(N_ELEMENTS(size_xbreak)-2) EQ  TYP_DOUBLE) THEN type = TYP_DOUBLE
   IF (KEYWORD_SET(double) EQ true) THEN type = TYP_DOUBLE
   ;
   ; Setup the parameters for the call to the system function.
   ;
   ; Floating point arguments and keywords   
   IF (type EQ TYP_DOUBLE) THEN BEGIN
      mach = imsl_machine(/DOUBLE)
      result1 = dblarr(npde, nx, n_t_vals)
      t_cvt = double(t)
      xbreak_cvt = double(xbreak)
      IF (KEYWORD_SET(tolerance)) THEN tolerance_cvt  =  DOUBLE(tolerance(0)) $
        ELSE tolerance_CVT = 100.*(mach.max_rel_space)
      IF (KEYWORD_SET(hinit)) THEN hinit_cvt  =  DOUBLE(hinit(0)) $
        ELSE hinit_CVT = DOUBLE(0.0)
   END ELSE BEGIN
      mach = imsl_machine(/FLOAT)
      result1 = FLTARR(npde, nx, n_t_vals)
      t_cvt = FLOAT(t)
      xbreak_cvt = FLOAT (xbreak)
      IF (KEYWORD_SET(tolerance)) THEN tolerance_cvt  =  FLOAT(tolerance(0)) $
        ELSE tolerance_CVT = 100.*(mach.max_rel_space)
      IF (KEYWORD_SET(hinit)) THEN hinit_cvt  =  FLOAT(hinit(0)) $
        ELSE hinit_CVT = FLOAT(0.0)
   END
   npde_cvt  =  IMSL_LONG(npde(0))
   IF (KEYWORD_SET(deriv_init)) THEN BEGIN
     deriv_init_cvt     =     imsl_cvt_arr(deriv_init,     type)
     IF (type EQ TYP_FLOAT) THEN  deriv_spc   =    FLTARR(npde,   nx,   n_t_vals) $
       ELSE  deriv_spc   =    DBLARR(npde,   nx,   n_t_vals)
      deriv_spc(*, *, 0)  =  deriv_init_cvt
   END
   y_cvt = imsl_cvt_arr(y, type)
   ;
   ; Call the system function.
   ;
   err_status = 0L
   MATHSTAT_224, type, err_status, npde_cvt, $
                              y_cvt, $
                              t_cvt, $
                              xbreak_cvt, $
                              nx, $
                              f_ut, $
                              f_bc, $
                              n_t_vals, $
                              hinit_cvt, $
                              tolerance_cvt, $
                              deriv_init_cvt, $
                              deriv_spc, $
                              result1

   ; Since the result of the routine is a 3D array, we have to
   ; adjust the order of the elements in the result to match
   ; the documented order.
   IF (type EQ TYP_DOUBLE) THEN result   =   DBLARR(npde_cvt,  nx,  n_t_vals) $
     ELSE  result   =   FLTARR(npde_cvt,  nx,  n_t_vals) 
   FOR i  =  0,   n_t_vals-1 DO BEGIN
      time_block   =   i*(npde_cvt*nx)
      FOR j = 0, npde_cvt-1 DO BEGIN
         result(j,  *,  i)   =   result1((time_block+j*nx):(time_block+(j+1)*nx-1))
      END
   END
   IF (KEYWORD_SET(deriv_init)) THEN BEGIN
      IF (type EQ TYP_DOUBLE) THEN deriv_init   =   DBLARR(npde_cvt,  nx,  n_t_vals) $
        ELSE  deriv_init   =   FLTARR(npde_cvt,  nx,  n_t_vals) 
      FOR i  =  0,   n_t_vals-1 DO BEGIN
         time_block   =   i*(npde_cvt*nx)
         FOR j = 0, npde_cvt-1 DO BEGIN
            deriv_init(j,  *,  i)   =   deriv_spc((time_block+j*nx):(time_block+(j+1)*nx-1))
         END
      END
   END
   RETURN,  result
END

