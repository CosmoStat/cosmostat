; $Id: //depot/idl/IDL_70/idl_src/libs/imsl/imsl_6.0/lib/imsl_radbf.pro#1 $
;
; Copyright (c) 1970-2006, VISUAL NUMERICS Inc. All Rights Reserved.
; This software is confidential information which is proprietary to and a
; trade secret of Visual Numerics Inc.  Use, duplication or disclosure is
; subject to the terms of an appropriate license agreement.
;
FUNCTION imsl_radbf, abscissa, $              ;INPUT 2-D array: floating point 
                fdata, $                   ;INPUT 1-D array: floating point 
                numcenters, $              ;INPUT Scalar LONG
                centers=centers, $         ;INPUT 2-D array: floating point 
                ratio_centers=ratio_centers, $ ;INPUT Scalar floating point
                random_seed=random_seed, $ ;INPUT Scalar LONG
                basis=basis, $             ;INPUT Scalar STRING
                delta=delta, $             ;INPUT Scalar floating point
                qr=qr, $                   ;INPUT Scalar floating point
                weights=weights, $         ;INPUT 1-D array: floating point
                double=double              ;INPUT Scalar ON/OFF flag

@imsl_init.pro
   ON_ERROR, on_err_action
   ;
   ; Error checking.
   ; - abscissa must be a 2-D array (size  ndim x ndata )
   ; - fdata must be a 1-D array (with ndata elements)
   ; - numcenters is converted to TYP_MEMINT, and must be positive.
   ; - If WEIGHTS set, it must be a 1-D array (with ndata elements)
   ; - If CENTERS set, it must be a 2-D array of size (numcenters x ndata)
   ; - If BASIS is supplied, then we will be using a user supplied WAVE
   ;   function, otherwise, we will use our default function.
   ;                       
   nargs = n_params()
   IF (nargs NE 3) THEN message, 'Incorrect number of arguments.'
   size_abs = IMSL_SIZE(abscissa)
   size_fdata = IMSL_SIZE(fdata)
   IF (size_abs(0) NE 2) THEN message, 'ABSCISSA must be a 2-D array.'
   ndim = IMSL_LONG(size_abs(1))
   ndata = IMSL_LONG(size_abs(2))
   IF ((size_fdata(0) NE 1) OR (N_ELEMENTS(fdata) NE ndata)) THEN $
     message, 'FDATA is not the correct size.'
   numcenters_cvt = IMSL_LONG(numcenters(0))
   IF (numcenters_cvt LT 1) THEN message, "The number of centers must be positive"
   IF (KEYWORD_SET(weights)) THEN BEGIN 
      size_weights = IMSL_SIZE(weights)
      IF ((size_weights(0) NE 1) OR (N_ELEMENTS(weights) NE ndata)) THEN $
        message, 'WEIGHTS must have the same dimensions as FDATA.'
   END
   IF (KEYWORD_SET(centers)) THEN BEGIN 
      size_centers = IMSL_SIZE(centers)
      IF ((size_centers(0) NE 1) OR (N_ELEMENTS(centers) NE numcenters_cvt)) THEN $
        message, 'CENTERS must be a 1-D array OF length NUMCENTERS.'
   END
   IF (KEYWORD_SET(weights)) THEN BEGIN 
      size_weights = IMSL_SIZE(weights)
      IF ((size_weights(0) NE 1) OR (N_ELEMENTS(weights) NE ndata)) THEN $
        message, 'WEIGHTS must have the same dimensions as XDATA.'
   END
   IF (KEYWORD_SET(basis)) THEN BEGIN
      size_basis = IMSL_SIZE(basis)
      IF ((N_ELEMENTS(basis) NE 1) OR (size_basis(N_ELEMENTS(size_basis)-2) NE 7)) THEN $
        message, 'BASIS must be a scalar string.'
   END
   ; 
   ; Decide on what precision to use.
   type = TYP_FLOAT
   IF (size_abs(N_ELEMENTS(size_abs)-2) EQ  TYP_DOUBLE) THEN type = TYP_DOUBLE
   IF (size_fdata(N_ELEMENTS(size_fdata)-2) EQ  TYP_DOUBLE) THEN type = TYP_DOUBLE
   IF (KEYWORD_SET(double) EQ true) THEN type = TYP_DOUBLE ELSE type = TYP_FLOAT
   ;
   ; Setup the parameters for the call to the system function.
   ; Input LONG keyword(s)
   ;
   IF (KEYWORD_SET(random_seed)) THEN random_seed_cvt = IMSL_LONG(random_seed(0)) $
     ELSE random_seed_cvt = IMSL_LONG(234579)
   IF (KEYWORD_SET(qr)) THEN qr_cvt = IMSL_1
   ;
   ; Floating point arguments and keywords   
   IF (type EQ TYP_DOUBLE) THEN BEGIN
      abscissa_cvt = double(abscissa)
      fdata_cvt = double(fdata)
      IF (KEYWORD_SET(centers)) THEN centers_cvt = double(centers)
      ; This bit of messy code is because 0.0 is a valid input value for RATIO_CENTERS.
      ratio_cvt = double(0.5)
      IF (arg_present(ratio_centers)) THEN $
        IF N_ELEMENTS(ratio_centers) THEN ratio_cvt = double(ratio_centers(0))
      IF (KEYWORD_SET(delta)) THEN delta_cvt = double(delta(0)) $
        ELSE delta_cvt = double(1.)
      IF (KEYWORD_SET(weights)) THEN weights_cvt = double(weights)
   END ELSE BEGIN
      abscissa_cvt = float(abscissa)
      fdata_cvt = float(fdata)
      IF (KEYWORD_SET(centers)) THEN centers_cvt = float(centers)
      ; This bit of messy code is because 0.0 is a valid input value for RATIO_CENTERS.
      ratio_cvt = float(0.5)
      IF (arg_present(ratio_centers)) THEN $
        IF N_ELEMENTS(ratio_centers) THEN ratio_cvt = float(ratio_centers(0))
      IF (KEYWORD_SET(delta)) THEN delta_cvt = float(delta(0)) $
        ELSE delta_cvt = float(1.)
      IF (KEYWORD_SET(weights)) THEN weights_cvt = float(weights)
   END
   ;
   ; Define the result variable. This variable will be filled with a
   ; radial basis structure if the computations end successfully.
   result = IMSL_1
   ;
   ; Call the system function.
   ;
   err_status = 0L
   MATHSTAT_180, type, err_status, abscissa_cvt, fdata_cvt, numcenters_cvt, $
                              ndim, ndata, $
                              centers_cvt, $
                              ratio_cvt, $
                              random_seed_cvt, $
                              basis, $
                              delta_cvt, $
                              weights_cvt, $
                              qr_cvt, $
                              result
 RETURN, result
END

