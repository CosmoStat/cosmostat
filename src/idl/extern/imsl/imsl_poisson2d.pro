; $Id: //depot/idl/IDL_70/idl_src/libs/imsl/imsl_6.0/lib/imsl_poisson2d.pro#1 $
;
; Copyright (c) 1970-2006, VISUAL NUMERICS Inc. All Rights Reserved.
; This software is confidential information which is proprietary to and a
; trade secret of Visual Numerics Inc.  Use, duplication or disclosure is
; subject to the terms of an appropriate license agreement.
;
FUNCTION imsl_poisson2d, rhs_pde, $            ;INPUT Scalar STRING
                    rhs_bc, $             ;INPUT Scalar STRING
                    coef_u, $             ;INPUT scalar floating point
                    nx, $                 ;INPUT Scalar long
                    ny, $                 ;INPUT Scalar long
                    ax, $                 ;INPUT Scalar floating point
                    bx, $                 ;INPUT Scalar floating point
                    ay, $                 ;INPUT Scalar floating point
                    by, $                 ;INPUT Scalar floating point
                    bc_type, $            ;INPUT 1-D array floating point
                    order=order, $        ;INPUT Scalar long.
                    double = DOUBLE       ;INPUT Scalar ON/OFF flag
                    
@imsl_init.pro
   ON_ERROR, on_err_action
   ;
   ; Error checking.
   ; - RHS_PDE must be a scalar string.
   ; - RHS_BC must be a scalar string.
   ; - BC_TYPE must be a 1-D array of length 4.
   ;                       
   nargs = n_params()
   IF (nargs NE 10) THEN message, 'Incorrect number of arguments.'
   size_rhs_pde = IMSL_SIZE(rhs_pde)
   IF ((N_ELEMENTS(rhs_pde) NE 1) OR (size_rhs_pde(N_ELEMENTS(size_rhs_pde)-2) NE 7)) THEN $
     message, 'RHS_PDE must be a scalar string.'
   size_rhs_bc = IMSL_SIZE(rhs_bc)
   IF ((N_ELEMENTS(rhs_bc) NE 1) OR (size_rhs_bc(N_ELEMENTS(size_rhs_bc)-2) NE 7)) THEN $
     message, 'RHS_BC must be a scalar string.'

   size_bc_type = IMSL_SIZE(bc_type)
      IF (size_bc_type(0) NE 1) THEN message, 'BC_TYPE must be a 1-D array.'
      IF (size_bc_type(1) NE 4) THEN message, 'BC_TYPE is not the correct size.'
   ; 
   ; Decide on what precision to use.
   type = TYP_FLOAT
   IF (KEYWORD_SET(DOUBLE) EQ true) THEN type  =  TYP_DOUBLE ELSE BEGIN 
      size_coef_u = IMSL_SIZE(coef_u)
      size_ax = IMSL_SIZE(ax)
      size_bx = IMSL_SIZE(bx)
      size_ay = IMSL_SIZE(ay)
      size_by = IMSL_SIZE(by)
      IF (size_coef_u(N_ELEMENTS(size_coef_u)-2) EQ  TYP_DOUBLE) THEN type  =  TYP_DOUBLE
      IF (size_ax(N_ELEMENTS(size_ax)-2) EQ  TYP_DOUBLE) THEN type  =  TYP_DOUBLE
      IF (size_bx(N_ELEMENTS(size_bx)-2) EQ  TYP_DOUBLE) THEN type  =  TYP_DOUBLE
      IF (size_ay(N_ELEMENTS(size_ay)-2) EQ  TYP_DOUBLE) THEN type  =  TYP_DOUBLE
      IF (size_by(N_ELEMENTS(size_by)-2) EQ  TYP_DOUBLE) THEN type   =   TYP_DOUBLE
   END
   ;
   ; Setup the parameters for the call to the system function.
   ; Input LONG keyword(s)
   ;
   IF (KEYWORD_SET(order)) THEN order_cvt = IMSL_LONG(order(0)) ELSE order_cvt = IMSL_4
   nx_cvt  =  IMSL_LONG(nx(0))
   ny_cvt  =  IMSL_LONG(ny(0))
   bc_type_cvt = IMSL_LONG(bc_type)
   ;
   ; Floating point arguments and keywords   
   IF (type EQ TYP_DOUBLE) THEN BEGIN
      result     = DBLARR(ny_cvt,nx_cvt)
      coef_u_cvt = DOUBLE(coef_u(0))
      ax_cvt = DOUBLE(ax(0))
      bx_cvt = DOUBLE(bx(0))
      ay_cvt = DOUBLE(ay(0))
      by_cvt  =DOUBLE(by(0))
   END ELSE BEGIN
      result     = FLTARR(ny_cvt,nx_cvt)
      coef_u_cvt = FLOAT(coef_u(0))
      ax_cvt = FLOAT(ax(0))
      bx_cvt = FLOAT(bx(0))
      ay_cvt = FLOAT(ay(0))
      by_cvt  =FLOAT(by(0))
   END
   ;
   ; Call the system function.
   ;
   err_status = 0L
   MATHSTAT_220, type, err_status, $
                         rhs_pde, $
                         rhs_bc, $
                         coef_u_cvt, $
                         ax_cvt, $
                         bx_cvt, $
                         ay_cvt, $
                         by_cvt, $
                         nx_cvt, $
                         ny_cvt, $
                         order_cvt, $
                         bc_type_cvt, $
                         result

 RETURN, TRANSPOSE(result)
END

