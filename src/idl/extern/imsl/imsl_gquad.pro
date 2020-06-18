; $Id: //depot/idl/IDL_70/idl_src/libs/imsl/imsl_6.0/lib/imsl_gquad.pro#1 $
;
; Copyright (c) 1970-2006, VISUAL NUMERICS Inc. All Rights Reserved.
; This software is confidential information which is proprietary to and a
; trade secret of Visual Numerics Inc.  Use, duplication or disclosure is
; subject to the terms of an appropriate license agreement.
;
PRO imsl_gquad, n, $                                 ;INPUT Scalar LONG
                  weights, $                      ;OUTPUT 1-D array: floating point
                  points, $                       ;OUTPUT 1-D array: floating point
                  Double=double, $                ;INPUT Scalar ON/OFF flag
                  cheby_first=cheby_first, $      ;INPUT Scalar ON/OFF flag
                  cheby_second=cheby_second, $    ;INPUT Scalar ON/OFF flag
                  cosh=cosh, $                    ;INPUT Scalar ON/OFF flag
                  fixed_points=fixed_points, $    ;INPUT 1-D array; floating point
                  hermite=hermite, $              ;INPUT Scalar ON/OFF flag
                  jacobi=jacobi, $                ;INPUT 1-D array; floating point
                  laguerre=laguerre               ;INPUT Scalar floating point
@imsl_init.pro
   ON_ERROR, on_err_action
   ;
   ; Error checking.
   ;   - n is converted to a long scalar. Must be GT 0.
   ;   - Only one weight choice is allowed.
   ;   - if JACOBI is set, it must specify an array of length 2.
   ;   - if FIXED_POINTS is set, it must specify an array of length 1 or 2.
   ;                           of a scalar, which is treated the same as 
   ;                           an array of length 1.
   ;
   nargs = n_params()
   IF (nargs NE 3) THEN $
     message, "Incorrect number of arguments."
   n_cvt = IMSL_LONG(n(0))
   IF (n_cvt LT 1) THEN message, 'N must be positive.'
   tmp = KEYWORD_SET(cosh) + KEYWORD_SET(fixed_points) + KEYWORD_SET(jacobi) + $
     KEYWORD_SET(laguerre) + KEYWORD_SET(cheby_first) + $
     KEYWORD_SET(cheby_second) + KEYWORD_SET(hermite)
   IF (tmp GT 1) THEN message, 'At most one of weight choice may be specified.'
   IF (KEYWORD_SET(jacobi)) THEN BEGIN
      size_jacobi = IMSL_SIZE(jacobi)
      IF ((size_jacobi(0) NE 1) OR (N_ELEMENTS(jacobi) NE 2)) THEN $
        message, 'JACOBI must be a 1-D array of length two.'
   END
   IF (KEYWORD_SET(fixed_points)) THEN BEGIN
      n_fixed_points = IMSL_N_ELEMENTS(fixed_points)
      IF (n_fixed_points  GT 2) THEN $
        message, 'FIXED_POINTS can have at most two elements.'
   END
   ;
   ; Decide on what precision to use.
   ; Use the highest precision of the input argument(s).
   ;
   type = TYP_FLOAT
   IF (KEYWORD_SET(double) EQ true) THEN type = TYP_DOUBLE
   ;
   ; Setup the parameters for the call to the system function.
   ;
   IF (type EQ TYP_DOUBLE) THEN BEGIN
      weights_spc = dblarr(n_cvt)
      points_spc = dblarr(n_cvt)
      IF (KEYWORD_SET(fixed_points) EQ TRUE) THEN fixed_pts_cvt = double(fixed_points)
      IF (KEYWORD_SET(jacobi) EQ TRUE) THEN jacobi_cvt = double(jacobi)
      IF (KEYWORD_SET(laguerre) EQ TRUE) THEN laguerre_cvt = double(laguerre)
   END ELSE BEGIN
      weights_spc = fltarr(n_cvt)
      points_spc = fltarr(n_cvt)
      IF (KEYWORD_SET(fixed_points) EQ TRUE) THEN fixed_pts_cvt = float(fixed_points)
      IF (KEYWORD_SET(jacobi) EQ TRUE) THEN jacobi_cvt = float(jacobi)
      IF (KEYWORD_SET(laguerre) EQ TRUE) THEN laguerre_cvt = float(laguerre)
   END
   IF (KEYWORD_SET(jacobi) EQ TRUE) THEN BEGIN
      jacobi_0 = jacobi_cvt(0)
      jacobi_1 = jacobi_cvt(1)
   END
   IF (KEYWORD_SET(fixed_points) EQ TRUE) THEN BEGIN 
      fixed_pts_0 = fixed_pts_cvt(0)
      IF (n_fixed_points EQ 2) THEN fixed_pts_1 = fixed_pts_cvt(1)
   END
   
   ;
   ; Call the system function.
   ;
   err_status = 0L
   MATHSTAT_148, type, err_status, n_cvt, weights_spc, points_spc, $
                  IMSL_LONG(KEYWORD_SET(cheby_first)), $
                  IMSL_LONG(KEYWORD_SET(cheby_second)), $
                  IMSL_LONG(KEYWORD_SET(cosh)), $
                  IMSL_LONG(KEYWORD_SET(hermite)), $
                  n_fixed_points, $
                  fixed_pts_0, $
                  fixed_pts_1, $
                  jacobi_0, $
                  jacobi_1, $
                  laguerre_cvt

   weights =  weights_spc
   points =  points_spc
   ;
   ; Return.
   ;
   RETURN
END
   

