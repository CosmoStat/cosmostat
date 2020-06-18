; $Id: //depot/idl/IDL_70/idl_src/libs/imsl/imsl_6.0/lib/imsl_convol1d.pro#1 $
;
; Copyright (c) 1970-2006, VISUAL NUMERICS Inc. All Rights Reserved.
; This software is confidential information which is proprietary to and a
; trade secret of Visual Numerics Inc.  Use, duplication or disclosure is
; subject to the terms of an appropriate license agreement.
;
FUNCTION imsl_Convol1d, x, $
                   y, $
                   Periodic=periodic, Direct=direct
;+
; NAME:
;       CONVOL1D
; PURPOSE:
;       Compute the discrete convolution of two arrays.
; CATEGORY:
;       Transforms, Signal Processing
; CALLING SEQUENCE:
;       Z = CONVOL1D(X, Y)
; INPUTS:
;       X   1-D array.
;       Y   1-D array.
; KEYWORD
;       Periodic:  If present and non-zero, then a periodic convolution
;                  is computed.
;       Direct:    If present and non-zero, then use the direct method
;                  (as opposed to the usually-quicker FFT) regardless
;                  of the size of the input.
; OUTPUTS:
;       1-D array containing the discrete convolution of X and Y.
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
   IF (nargs NE 2) THEN  message, 'Incorrect number of arguments.'
   ;
   ; Error checking:
   ;    X must be a 1-D array
   ;    Y must be a 1-D array
   size_x = IMSL_SIZE(x)
   IF (size_x(0) NE 1) THEN message, 'X must be a 1-D array.'
   nx = IMSL_N_ELEMENTS(x)
   size_y = IMSL_SIZE(y)
   IF (size_y(0) NE 1) THEN message, 'Y must be a 1-D array.'
   ny = IMSL_N_ELEMENTS(y)
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
   nz = imsl_POW235(nx, ny, ipad)
   ;
   ; Compute length of zero padding, based upon NZ.
   nxpad = nz-nx
   nypad = nz-ny
   ;
   ; Convert the data, and zero-pad at the end of the vectors.
   ; We need to convert it since we will use a slightly different algorithm
   ; for complex data than real data.
   IF (type LT TYP_COMPLEX) THEN BEGIN
      IF (nxpad GT 0) THEN x_cvt = [x, dblarr(nxpad)] ELSE x_cvt = double(x)
      IF (nypad GT 0) THEN y_cvt = [y, dblarr(nypad)] ELSE y_cvt = double(y)
      zh = dblarr(nz)  ; Always return double precision.
   END ELSE BEGIN
      IF (nxpad GT 0) THEN x_cvt = [x, dcomplexarr(nxpad)] ELSE x_cvt = x
      IF (nypad GT 0) THEN y_cvt = [y, dcomplexarr(nypad)] ELSE y_cvt = y
   END
   ;
   ; Now x_cvt and y_cvt are the correct data type and are padded correctly.
   ;
   ; If IPAD is 1, and the length of the result is less than NDIRECT,
   ; then compute the convolution directly.
   ndirect = 16L
   IF (((ipad EQ 1) AND (nz LT ndirect)) OR KEYWORD_SET(direct)) THEN BEGIN
      z = make_array(nx + ny -1, type = IMSL_SIZE(x_cvt,/type))
      FOR i = IMSL_0, nx-1 DO BEGIN
         FOR j = IMSL_0, ny-1 DO BEGIN
            k = i+j
            z(k) = z(k) + x(i)*y(j)
         ENDFOR
      ENDFOR
      RETURN, z
   END
   ; In the case that x and y are real, we compute the ffts using the
   ; real transform in double precision.
   IF (type NE TYP_COMPLEX AND type NE TYP_DCMPLX) THEN BEGIN
      ; Real forward transform the padded X and Y.
      xh = imsl_fftcomp(x_cvt, /DOUBLE)
      yh = imsl_fftcomp(y_cvt, /DOUBLE)
      ; Compute pre-reverse-transform real vector. This mimicks
      ; the process of taking the complex forward transform of the
      ; original vectors, multiplying the elements of the transformed
      ; vectors, and converting that back to the real vector to be
      ; reverse-transformed.
      zh(0) = xh(0)*yh(0)
      idx = lindgen((nz-1)/2)
      idx1 = 2*idx+1
      idx2 = 2*(idx+1)
      zh(idx1) = xh(idx1)*yh(idx1) - xh(idx2)*yh(idx2)
      zh(idx2) = xh(idx2)*yh(idx1) + xh(idx1)*yh(idx2)
      IF ((nz MOD 2) EQ 0) THEN zh(nz-1) = xh(nz-1)*yh(nz-1)
      ; Reverse transform real the vector and scale the
      ; answer by NZ.
      RETURN, imsl_fftcomp(zh, /backward, /DOUBLE)/nz
   END ELSE BEGIN
      ; In the complex case, we want to have precision
      ; equivalent to double precision, but WAVE only supports single
      ; precision complex.  So, we use the routine DCMPLXFFT, which
      ; uses the double precision real FFT to compute a double
      ; precision complex FFT stored in seperate real and imaginary
      ; arrays.
      ; Transform X.
      imsl_dcmplxfft, double(x_cvt), imaginary(x_cvt), xhr, xhi
      imsl_dcmplxfft, double(y_cvt), imaginary(y_cvt), yhr, yhi
      ;
      ; Multiply component-wise.  Note were are multiplying
      ; complex numbers.
      zhr = xhr*yhr - xhi*yhi
      zhi = xhr*yhi + xhi*yhr
      ;
      ; Compute inverse FFT of complex(zhr, zhi) in double precision
      ; using DCMPLXFFT().
      imsl_dcmplxfft, zhr, zhi, result_r, result_i, /BACKWARD
      ;
      IF type EQ TYP_COMPLEX THEN RETURN, COMPLEX(result_r,result_i)/nz ELSE $
                                  RETURN, DCOMPLEX(result_r,result_i)/nz
   END
END
