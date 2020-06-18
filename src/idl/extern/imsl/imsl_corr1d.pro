; $Id: //depot/idl/IDL_70/idl_src/libs/imsl/imsl_6.0/lib/imsl_corr1d.pro#1 $
;
; Copyright (c) 1970-2006, VISUAL NUMERICS Inc. All Rights Reserved.
; This software is confidential information which is proprietary to and a
; trade secret of Visual Numerics Inc.  Use, duplication or disclosure is
; subject to the terms of an appropriate license agreement.
;
FUNCTION imsl_Corr1d, x, $
                 y, $
                 PERIODIC=periodic
;+
; NAME:
;       CORR1D
; PURPOSE:
;       Compute the discrete correlation of two 1-D arrays.
; CATEGORY:
;       Transforms, Signal Processing
; CALLING SEQUENCE:
;       Z = CORR1D(X[, Y])
; INPUTS:
;       X   1-D array.
;       Y   1-D array.
; KEYWORD
;       PERIODIC:  If present and non-zero, then a periodic correlation
;                  is computed.
; OUTPUTS:
;       1-D array containing the discrete correlation of X and Y.
;       If only one argument is supplied then the discrete
;       correlationof X and X is returned.
; MODIFICATION HISTORY:
;       6, May 1994, Written by Mike Pulverenti for VNI.
;-
@imsl_init.pro
   ON_ERROR, on_err_action
   ;
   ; Value of IPAD is used to determine amount of zero-padding.
   IF KEYWORD_SET(periodic) THEN ipad = 0 ELSE ipad = 1
   ;
   nargs = n_params()
   IF ((nargs NE 2) AND (nargs NE 1)) THEN  $
     message, 'Incorrect number of arguments.'
   ;
   ; Error checking:
   ;    X must be a 1-D array
   ;    Y must be a 1-D array, if present, and must be same length as X.
   size_x = IMSL_SIZE(x)
   IF (size_x(0) NE 1) THEN message, 'X must be a 1-D array.'
   nx = IMSL_N_ELEMENTS(x)
   IF (nargs EQ 2) THEN BEGIN
      size_y = IMSL_SIZE(y)
      IF (size_y(0) NE 1) THEN message, 'Y must be a 1-D array.'
      ny = IMSL_N_ELEMENTS(y)
      IF (ny NE nx) THEN message, 'Y must be the same length as X.'
   END ELSE BEGIN
      ; These values are defined for convenience.
      size_y = size_x
      ny = nx
   END
   ;
   ; Decide on what precision to use.
   ; Decision based upon highest type of X and Y.  Lowest precision
   ; used is TYP_FLOAT, highest precision used is TYP_DCMPLX.
   ;
   type = TYP_FLOAT
   highest_type = size_x(N_ELEMENTS(size_x)-2) > size_y(N_ELEMENTS(size_y)-2)
   IF (highest_type GT TYP_FLOAT) THEN type = TYP_DOUBLE
   IF (highest_type GT TYP_DOUBLE) THEN type = TYP_COMPLEX
   IF (highest_type GT TYP_COMPLEX) THEN type = TYP_DCMPLX
   ;
   ; Compute the length of the returned vector Z.
   nz = IMSL_POW235(nx, ny, ipad)
   ;
   ; Compute length of zero padding, based upon NZ.
   nxpad = nz-nx
   nypad = nz-ny
   ;
   ; Convert the data, and zero-pad at the end of the vectors.
   ; We need to convert it since we will use a slightly different algorithm
   ; for complex data than real data.
   IF (type LT TYP_COMPLEX) THEN BEGIN
      IF (nxpad GT 0) THEN x_cvt = [x, fltarr(nxpad)] ELSE x_cvt = x
      IF (nargs EQ 2) THEN BEGIN
         IF (nypad GT 0) THEN y_cvt = [y, fltarr(nypad)] ELSE y_cvt = y
      END
      IF (type EQ TYP_DOUBLE) THEN zh = dblarr(nz) ELSE zh = fltarr(nz)
   END ELSE BEGIN
      IF (nxpad GT 0) THEN x_cvt = [x, dcomplexarr(nxpad)] ELSE x_cvt = x
      IF (nargs EQ 2) THEN begin
         IF (nypad GT 0) THEN y_cvt = [y, dcomplexarr(nypad)] ELSE y_cvt = y
      END

   END
   ;
   ; Now x_cvt and y_cvt are the correct data type and are padded correctly.
   ;
   ; Seperate cased for one and two positional arguments.
   ;
   IF (nargs EQ 2) THEN BEGIN
      IF (type LT TYP_COMPLEX) THEN BEGIN
         ; Real forward transform the padded X and Y.
         xh = imsl_fftcomp(x_cvt)
         yh = imsl_fftcomp(y_cvt)
         ; Compute pre-reverse-transform real vector. This mimicks
         ; the process of taking the complex forward transform of the
         ; original vectors, multiplying the conjugates elements of the
         ; transformed vectors, and converting that back to the real vector
         ; to be reverse-transformed.
         zh(0) = xh(0)*yh(0)
         idx = lindgen((nz-1)/2)
         idx1 = 2*idx+1
         idx2 = 2*(idx+1)
         zh(idx1) = xh(idx1)*yh(idx1) + xh(idx2)*yh(idx2)
         zh(idx2) = xh(idx2)*yh(idx1) - xh(idx1)*yh(idx2)
         IF ((nz MOD 2) EQ 0) THEN zh(nz-1) = xh(nz-1)*yh(nz-1)
         ; Reverse transform real the vector and scale the
         ; answer by NZ.
         RETURN, imsl_fftcomp(zh, /backward)/nz
      END ELSE BEGIN
         ; Complex forward transform the padded X and Y, then
         ; backward transform and scale the result of the element-wise
         ; multiplication of the transformed X and CONJ(Y).
         RETURN, imsl_fftcomp(imsl_fftcomp(x_cvt)*conj(imsl_fftcomp(y_cvt)), /back)/nz
      END
   END ELSE BEGIN
      IF (type LT TYP_COMPLEX) THEN BEGIN
         ; Real forward transform the padded X.
         xh = imsl_fftcomp(x_cvt)
         ; Compute pre-reverse-transform real vector. This mimicks
         ; the process of taking the complex forward transform of the
         ; original vectors, multiplying the conjugates elements of the
         ; transformed vectors, and converting that back to the real vector
         ; to be reverse-transformed.
         zh(0) = xh(0)*xh(0)
         idx = lindgen((nz-1)/2)
         idx1 = 2*idx+1
         idx2 = 2*(idx+1)
         zh(idx1) = xh(idx1)*xh(idx1) + xh(idx2)*xh(idx2)
         zh(idx2) = xh(idx2)*xh(idx1) - xh(idx1)*xh(idx2)
         IF ((nz MOD 2) EQ 0) THEN zh(nz-1) = xh(nz-1)*xh(nz-1)
         ; Reverse transform real the vector and scale the
         ; answer by NZ.
         RETURN, imsl_fftcomp(zh, /backward)/nz
      END ELSE BEGIN
         ; Complex forward transform the padded X, then
         ; backward transform and scale the result of the element-wise
         ; multiplication of the transformed X and CONJ(transformed(X)).
         RETURN, imsl_fftcomp(imsl_fftcomp(x_cvt)*conj(imsl_fftcomp(x_cvt)), /back)/nz
      END
   END
END
