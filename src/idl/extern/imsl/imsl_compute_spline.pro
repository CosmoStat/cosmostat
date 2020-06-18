; $Id: //depot/idl/IDL_70/idl_src/libs/imsl/imsl_6.0/lib/imsl_compute_spline.pro#1 $   
;
; Copyright (c) 1970-2006, VISUAL NUMERICS Inc. All Rights Reserved.
; This software is confidential information which is proprietary to and a
; trade secret of Visual Numerics Inc.  Use, duplication or disclosure is
; subject to the terms of an appropriate license agreement.
;
FUNCTION imsl_compute_spline, data1, $          ;INPUT 1-D array: floating point 
                data2, $                   ;INPUT 1-D array: floating point 
                data3, $                   ;INPUT 2-D array: floating point 
                fcn_idx=fcn_idx, $         ;INPUT Scalar LONG
                xorder=xorder, $           ;INPUT Scalar LONG
                yorder=yorder, $           ;INPUT Scalar LONG
                double=double, $           ;INPUT Scalar ON/OFF flag
                xknots=xknots, $           ;INPUT 1-D array: floating point
                yknots=yknots, $           ;INPUT 1-D array: floating point
                fcn_name=fcn_name, $       ;INPUT Scalar STRING.
                bspline=bspline, $         ;INPUT Scalar ON/OFF flag
                iright=iright, $           ;INPUT Scalar LONG
                right=right, $             ;INPUT Scalar floating point
                ileft=ileft, $             ;INPUT Scalar LONG
                left=left, $               ;INPUT Scalar floating point
                periodic=periodic, $       ;INPUT Scalar ON/OFF flag
                concave=concave, $         ;INPUT Scalar ON/OFF flag
                itmax=itmax, $             ;INPUT Scalar LONG
                sse=sse, $                 ;OUTPUT Scalar FLOAT
                xweights=xweights, $       ;INPUT 1-D array: floating point
                yweights=yweights, $       ;INPUT 1-D array: floating point
                optimize=optimize, $       ;INPUT Scalar ON/OFF flag
                xspace_dim=xspace_dim, $   ;INPUT Scalar LONG
                yspace_dim=yspace_dim, $   ;INPUT Scalar LONG
                smpar=smpar, $             ;INPUT Scalar floating point
                weights=weights            ;INPUT 1-D array: floating point 

@imsl_init.pro
   ON_ERROR, on_err_action
   ;
   ; Error checking:
   ;  Two cases: nargs = 2 or nargs = 3
   ;
   ;  CASE 1: nargs = 2: Compute 1-D spline.
   ;    - DATA1 must be 1-D. Set nx = N_ELEMENTS(DATA1)
   ;    - DATA2 must be 1-D array of length (nx). 
   ;    - YKNOTS, YORDER, YWEIGHTS not allowed.
   ;    - IRIGHT and RIGHT must be supplied together.
   ;    - ILEFT  and LEFT  must be supplied together.
   ;    - ITMAX is valid only if CONCAVE is also set.
   ;    - If WEIGHTS is present, it must be 1-D array of length (nx). 
   ;    The following checks are made after the CASE 2 so that
   ;    they only have to be coded once.
   ;    - If XORDER supplied, it must be positive.
   ;    - If XKNOTS supplied, it must be a 1-D array of length (nx + xorder).
   ;    - If XWEIGHTS are supplied, it must be a 1-D array of length (nx)
   ;
   ;  CASE 2: nargs = 3: Compute 2-D spline.
   ;    - DATA1 must be 1-D. Set nx = N_ELEMENTS(DATA1)
   ;    - DATA2 must be 1-D. Set ny = N_ELEMENTS(DATA2)
   ;    - DATA3 must be 2-D array of size (nx X ny)
   ;    - if XORDER supplied, it must be positive.
   ;    - if YORDER supplied, it must be positive.
   ;    - XKNOTS and YKNOTS must be supplied as a pair.
   ;    - XWEIGHTS and YWEIGHTS must be supplied as a pair.
   ;    - if XKNOTS supplied, it must be a 1-D array of length (nx + xorder). 
   ;    - if YKNOTS supplied, it must be a 1-D array of length (ny + xorder). 
   ;    - If XWEIGHTS are supplied, it must be a 1-D array of length (nx)
   ;    - If YWEIGHTS are supplied, it must be a 1-D array of length (ny)
   ;    - The keyword OPTIMIZE is not allowed in 2-D case.
   ;
   nargs = n_params()
   IF ((nargs NE 2) AND (nargs NE 3)) THEN message, 'Incorrect number of arguments.'
   size_data1 = IMSL_SIZE(data1)
   size_data2 = IMSL_SIZE(data2)
   IF (size_data1(0) NE 1) THEN message, 'XDATA must be a 1-D array.' $
     ELSE nx = IMSL_N_ELEMENTS(data1)
   
   IF (nargs EQ 2) THEN BEGIN
      IF ((size_data2(0) NE 1) OR (N_ELEMENTS(data2) NE nx)) THEN $
        message, 'FDATA must have the same dimensions as XDATA.'
      IF (KEYWORD_SET(yorder) + KEYWORD_SET(yknots) + KEYWORD_SET(yweights)) THEN $
        message, 'The keywords YORDER, YKNOTS or YWEIGHTS are not valid when ' + $
        'computing a 1-D spline.'
      domain_dim = IMSL_1
      IF (KEYWORD_SET(itmax)) THEN $
        IF (NOT KEYWORD_SET(concave)) THEN $
        message, "ITMAX is valid only IF CONCAVE is also set."
      IF (KEYWORD_SET(weights)) THEN BEGIN 
         size_weights = IMSL_SIZE(weights)
         IF ((size_weights(0) NE 1) OR (N_ELEMENTS(weights) NE nx)) THEN $
           message, 'WEIGHTS must have the same dimensions as XDATA.'
      END

      nycoefs    = IMSL_0
      yorder_cvt = IMSL_0
      ny         = IMSL_0
      yorder_cvt = IMSL_0
   END ELSE BEGIN
      ; (nargs eq 3)
      domain_dim = IMSL_2
      size_data3 = IMSL_SIZE(data3)
      IF (size_data2(0) NE 1) THEN message, 'YDATA must be a 1-D array.' $
        ELSE ny = IMSL_N_ELEMENTS(data2)
      IF (size_data3(0) NE 2) THEN $
        message, 'FDATA must be a 2-D array when computing a 2-D spline.'
      IF ((size_data3(1) NE nx) OR (size_data3(2) NE ny)) THEN $
        message, 'The dimensions of FDATA do not match those of XDATA nd YDATA.'
      IF (KEYWORD_SET(yorder)) THEN yorder_cvt = IMSL_LONG(yorder(0)) ELSE yorder_cvt = IMSL_4
      IF (yorder_cvt LT 1) THEN message, 'YORDER must be positive'
      IF (KEYWORD_SET(yspace_dim)) THEN nycoefs = yspace_dim ELSE nycoefs = ny
      IF (KEYWORD_SET(yknots)) THEN BEGIN
         IF (NOT (KEYWORD_SET(xknots))) THEN $
           message, 'The keywords XKNOTS and YKNOTS must be supplied together.'
         size_yknots = IMSL_SIZE(yknots)
         IF ((size_yknots(0) NE 1) OR (N_ELEMENTS(yknots) NE nycoefs+yorder_cvt)) THEN $
           message, 'YKNOTS must be a 1-D array of length ' + $
           strtrim(nycoef + yorder_cvt, 1) + '.'
      END
        
      IF (KEYWORD_SET(yweights)) THEN BEGIN
         IF (NOT (KEYWORD_SET(xweights))) THEN $
           message, 'The keywords XWEIGHTS and YWEIGHTS must be supplied together.'
         size_yweights = IMSL_SIZE(yweights)
         IF ((size_weights(0) NE 1) OR (N_ELEMENTS(yweights) NE ny)) THEN $
           message, 'YWEIGHTS must be a 1-D array of length ' + strtrim(ny, 1) + '.'
      END
      IF (KEYWORD_SET(optimize)) THEN $
        message, "The keyword OPTIMIZE is not allowed when computing a 2-D spline."
   END
   ;
   ; These checks are independent upon nargs.
   xorder_cvt = IMSL_4
   IF (KEYWORD_SET(xorder)) THEN xorder_cvt = IMSL_LONG(xorder(0))
   IF (xorder_cvt LT 1) THEN message, 'XORDER must be positive'
   IF (KEYWORD_SET(xspace_dim)) THEN nxcoefs = xspace_dim ELSE nxcoefs = nx
   IF (KEYWORD_SET(xknots)) THEN BEGIN
      size_xknots = IMSL_SIZE(xknots)
      IF ((size_xknots(0) NE 1) OR (N_ELEMENTS(xknots) NE nxcoefs+xorder_cvt)) THEN $
        message, 'XKNOTS must be a 1-D array of length ' + $
        strtrim(nxcoefs + xorder_cvt, 1) + '.'
   END
   IF (KEYWORD_SET(xweights)) THEN BEGIN
      size_xweights = IMSL_SIZE(xweights)
      IF ((size_weights(0) NE 1) OR (N_ELEMENTS(xweights) NE nx)) THEN $
        message, 'XWEIGHTS must be a 1-D array of length ' + strtrim(nx, 1) + '.'
   END
   total_knots = (nxcoefs + xorder_cvt) + (nycoefs + yorder_cvt) 
   total_coefs = nxcoefs * (nycoefs > IMSL_1) ;
   target_dim = IMSL_1
   
   ; Error checking complete.
   ;
   ; Decide on what precision to use.
   ;
   type = TYP_FLOAT
   IF (size_data1(N_ELEMENTS(size_data1)-2) EQ  TYP_DOUBLE) THEN type = TYP_DOUBLE
   IF (size_data2(N_ELEMENTS(size_data2)-2) EQ  TYP_DOUBLE) THEN type = TYP_DOUBLE
   IF (nargs EQ 3) THEN $
     IF (size_data3(N_ELEMENTS(size_data3)-2) EQ  TYP_DOUBLE) THEN type = TYP_DOUBLE
   IF (KEYWORD_SET(double) EQ true) THEN type = TYP_DOUBLE
   ;
   ; Setup the parameters for the call to the system function.
   ;----------------------------------------------------------
   ; Input LONG keyword(s)
   ;

   IF (KEYWORD_SET(bspline))  THEN compute_bspline = IMSL_1
   IF (KEYWORD_SET(iright)) THEN iright_cvt = IMSL_LONG(iright(0))
   IF (KEYWORD_SET(ileft))  THEN ileft_cvt  = IMSL_LONG(ileft(0))
   IF (KEYWORD_SET(periodic))  THEN periodic_cvt  = IMSL_LONG(periodic(0))
   IF (KEYWORD_SET(concave)) THEN concave_cvt  = IMSL_1
   IF (KEYWORD_SET(itmax)) THEN itmax_cvt  = IMSL_LONG(itmax(0))
   ; Output LONG keyword(s)
   ;
   ; Floating point arguments and keywords
   IF (type EQ TYP_DOUBLE) THEN BEGIN
      ;
      ; Input arguments and keyword(s)
      ;
      data1_cvt = double(data1)
      data2_cvt = double(data2)
      IF (nargs EQ 3) THEN data3_cvt = double(transpose(data3))
      IF (KEYWORD_SET(xknots)) THEN xknots_cvt = double(xknots)
      IF (KEYWORD_SET(yknots)) THEN yknots_cvt = double(yknots)
      IF (KEYWORD_SET(xweights)) THEN xweights_cvt = double(xweights)
      IF (KEYWORD_SET(yweights)) THEN yweights_cvt = double(yweights)
      ; Test iright and ileft since 0.0 is a valid value of right and left.
      IF (KEYWORD_SET(iright))  THEN right_cvt  = double(right(0))
      IF (KEYWORD_SET(ileft)) THEN left_cvt  = double(left(0))
      IF (KEYWORD_SET(weights)) THEN weights_cvt = double(weights)
      IF (KEYWORD_SET(smpar)) THEN smpar_cvt = double(smpar(0))
      ; Output keyword(s)
      ;
      IF (ARG_PRESENT(sse)) THEN sse_spc = double(0.0)
   END ELSE BEGIN
      ; Input arguments and keyword(s)
      ;
      data1_cvt = float(data1)
      data2_cvt = float(data2)
      IF (nargs EQ 3) THEN data3_cvt = float(transpose(data3))
      IF (KEYWORD_SET(xknots)) THEN xknots_cvt = float(xknots)
      IF (KEYWORD_SET(yknots)) THEN yknots_cvt = float(yknots)
      IF (KEYWORD_SET(xweights)) THEN xweights_cvt = float(xweights)
      IF (KEYWORD_SET(yweights)) THEN yweights_cvt = float(yweights)
      ; Test iright and ileft since 0.0 is a valid value of right and left.
      IF (KEYWORD_SET(iright))  THEN right_cvt  = float(right(0))
      IF (KEYWORD_SET(ileft))  THEN left_cvt  = float(left(0))
      IF (KEYWORD_SET(weights)) THEN weights_cvt = float(weights)
      IF (KEYWORD_SET(smpar)) THEN smpar_cvt = float(smpar(0))
      ; Output keyword(s)
      ;
      IF (ARG_PRESENT(sse)) THEN sse_spc = float(0.0)
   END
   ;
   ; Define the result variable. This variable will be filled with a
   ; spline structure if the computations end successfully.
; PPOLY
;   DOMAIN_DIM      LONG                 1
;   TARGET_DIM      LONG                 1
;   ORDER           LONG      Array(1)
;   NUM_COEF        LONG      Array(1)
;   NUM_BREAKPOINTS LONG      Array(1)
;   BREAKPOINTS     FLOAT     Array(10)
;   COEF            FLOAT     Array(36)
;BSPLINE
;   DOMAIN_DIM LONG                 1
;   TARGET_DIM LONG                 1
;   ORDER      LONG      Array(1)
;   NUM_COEF   LONG      Array(1)
;   NUM_KNOTS  LONG      Array(1)
;   KNOTS      FLOAT     Array(9)
;   COEF       FLOAT     Array(5)

   result1 = IMSL_1 ; DOMAIN_DIM
   result2 = IMSL_1 ; TARGET_DIM
   result3 = IMSL_1 ; ORDER
   result4 = IMSL_1 ; NUM_COEF
   result5 = IMSL_1 ; NUM_BREAKPOINTS / NUM_KNOTS
   result6 = IMSL_1 ; KNOTS/BREAKPOINTS
   result7 = IMSL_1 ; COEF
   ;
   ; Call the system function.
   ;
   err_status = 0L
   MATHSTAT_113, type, err_status, data1_cvt, data2_cvt, data3_cvt, nx, ny, $
                              xorder_cvt, yorder_cvt, $
                              xknots_cvt, yknots_cvt, $
                              sse_spc, $
                              xweights_cvt, $
                              yweights_cvt, $
                              xspace_dim, $
                              yspace_dim, $
                              optimize, $
                              iright_cvt, $
                              right_cvt, $
                              ileft_cvt, $
                              left_cvt, $
                              periodic_cvt, $
                              concave_cvt, $
                              itmax_cvt, $
                              weights_cvt, $
                              smpar_cvt, $
                              domain_dim, $
                              target_dim, $
                              fcn_name, $
                              fcn_idx-IMSL_1, $
                              compute_bspline, $
                              result1, $
                              result2, $
                              result3, $
                              result4, $
                              result5, $
                              result6, $
                              result7

   IF (ARG_PRESENT(sse)) THEN sse = sse_spc

   ; return
   if (n_elements(compute_bspline)) then begin
	result = { DOMAIN_DIM:result1, TARGET_DIM:result2, ORDER:result3, NUM_COEF:result4, NUM_KNOTS:result5, KNOTS:result6, COEF:result7}
   end else begin
	result = { DOMAIN_DIM:result1, TARGET_DIM:result2, ORDER:result3, NUM_COEF:result4, NUM_BREAKPOINTS:result5, BREAKPOINTS:result6, COEF:result7}
   end

   RETURN, result
END

                   
                   
                   

  
      

  
