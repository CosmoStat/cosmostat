; $Id: //depot/idl/IDL_70/idl_src/libs/imsl/imsl_6.0/lib/imsl_conlsq.pro#1 $   
;
; Copyright (c) 1970-2006, VISUAL NUMERICS Inc. All Rights Reserved.
; This software is confidential information which is proprietary to and a
; trade secret of Visual Numerics Inc.  Use, duplication or disclosure is
; subject to the terms of an appropriate license agreement.
;
function imsl_Chk_cnst_struct, constraints, type, num_constr
@imsl_init.pro
   ; The following checks are made on the argument that
   ; is supposed to be an array of constraint structures.
   ; - Must be an array of structures.
   ; - Must have 5 tags
   ; - The 1-st tag determines if the precision:TYP_FLOAT or TYP_DOUBLE
   ; - Tagnames must match expected tagnames.
   ; - For each tag:
   ;       o  check the data type.
   ;       o  check the size.
   ;
   err_str = "CONSTRAINTS must be a valid array of constraint structures."
   size_constr = IMSL_SIZE(constraints)
   IF ((size_constr(0) NE 1) OR $
       (size_constr(N_ELEMENTS(size_constr)-2) NE 8)) THEN $
     message, err_str
 
   n_tags = n_tags(constraints)
   IF (n_tags NE 5) THEN message, err_str
   constr_tagnames = tag_names(constraints)
   exp_tagnames = ["XVAL", "DER", "TYPE", "BL", "BU"]
   type = SIZE(constraints.(0), /TYPE)
   val = FIX(0, TYPE=type)
   num_constr = IMSL_N_ELEMENTS(constraints)
   struct = {XVAL: val, DER: IMSL_0, TYPE: IMSL_0, BL: val, BU: val}
   constraints_cvt = REPLICATE(struct, num_constr)
   for i = 0, n_tags-1 do begin
      constraints_cvt.(i) = constraints.(i)
   endfor

   RETURN, constraints_cvt
END

   
FUNCTION imsl_conlsq, xdata, $                ;INPUT 1-D array: floating point 
                fdata, $                   ;INPUT 1-D array: floating point 
                spacedim, $                ;INPUT Scalar LONG
                constraints, $             ;INPUT Array of structures
                nhard, $                   ;INPUT Scalar LONG
                weights=weights, $         ;INPUT 1-D array: floating point
                order=order, $             ;INPUT Scalar LONG
                double=double, $           ;INPUT Scalar ON/OFF flag
                knots=knots                ;INPUT 1-D array: floating point

@imsl_init.pro
   ;
   ; Error checking.
   ;    - xdata must be 1-D. Set nx = n_elements(xdata)
   ;    - fdata must be 1-D array of length (nx).
   ;    - spacedim is converted to a TYP_MEMINT scalar.
   ;    - constraints must be an array of stuctures.  The structures within the
   ;      array are tested seperately.
   ;    - if nargs == 5, nhard is converted to a TYP_MEMINT scalar.
   ;      array are checked to make sure they conform to what is expected.
   ;    - if ORDER supplied, it must be positive.
   ;    - if KNOTS supplied, it must be a 1-D array of length (space_dim + order). 
   ;    - if WEIGHTS supplied, it must be a 1-D array of length (nx). 
   ;                       
   ON_ERROR, on_err_action
   nargs = n_params()
   IF ((nargs NE 4) AND (nargs NE 5)) THEN message, 'Incorrect number of arguments.'
   size_xdata = IMSL_SIZE(xdata)
   size_fdata = IMSL_SIZE(fdata)
   IF (size_xdata(0) NE 1) THEN message, 'XDATA must be a 1-D array.' $
     ELSE nx = IMSL_N_ELEMENTS(xdata)
   IF ((size_fdata(0) NE 1) OR (N_ELEMENTS(fdata) NE nx)) THEN $
     message, 'FDATA must have the same dimensions as XDATA.'
   spacedim_cvt = IMSL_LONG(spacedim)
   IF (nargs EQ 5) THEN nhard_cvt = IMSL_LONG(nhard) ELSE nhard_cvt = IMSL_0
   ;
   ; Test out the constraints array.
   ;
   constraints_cvt = imsl_chk_cnst_struct(constraints, type, num_constr)

   ; CALL A FUNCTION TO TEST THE CONSTRAINTS ARRAY!!!!
   ;
   order_cvt = IMSL_4
   IF (KEYWORD_SET(order)) THEN order_cvt = IMSL_LONG(order(0))
   IF (order_cvt LT 1) THEN message, 'ORDER must be positive'
   IF (KEYWORD_SET(knots)) THEN BEGIN
      size_knots = IMSL_SIZE(knots)
      IF ((size_knots(0) NE 1) OR (N_ELEMENTS(knots) NE spacedim_cvt+order_cvt)) THEN $
        message, 'KNOTS must be a 1-D array of length ' + $
        strtrim(spacedim_cvt + order_cvt, 1) + '.'
   END
   IF (KEYWORD_SET(weights)) THEN BEGIN 
      size_weights = IMSL_SIZE(weights)
      IF ((size_weights(0) NE 1) OR (N_ELEMENTS(weights) NE nx)) THEN $
        message, 'WEIGHTS must have the same dimensions as XDATA.'
   END
   ; Decide on what precision to use.
   ; The precision used is determined by the input constraint array.
   ;
   ; Setup the parameters for the call to the system function.
   ; Floating point arguments and keywords
   IF (type EQ TYP_DOUBLE) THEN BEGIN
      ;
      ; Input arguments and keyword(s)
      ;
      xdata_cvt = double(xdata)
      fdata_cvt = double(fdata)
      IF (KEYWORD_SET(knots)) THEN knots_cvt = double(knots)
      IF (KEYWORD_SET(weights)) THEN weights_cvt = double(weights)
   END ELSE BEGIN
      ;
      ; Input arguments and keyword(s)
      ;
      xdata_cvt = float(xdata)
      fdata_cvt = float(fdata)
      IF (KEYWORD_SET(knots)) THEN knots_cvt = float(knots)
      IF (KEYWORD_SET(weights)) THEN weights_cvt = float(weights)
   END
   ;
   ; Define the result variables. These variable will be used to fill
   ; a spline structure if the computations end successfully.
   result1 = IMSL_1
   result2 = IMSL_1
   result3 = IMSL_1
   result4 = IMSL_1
   result5 = IMSL_1
   result6 = IMSL_1
   result7 = IMSL_1
   ;
   ; Call the system function.
   ;
   err_status = 0L
   MATHSTAT_123, type, err_status, xdata_cvt, fdata_cvt, nx, $
                              spacedim_cvt,$
                              constraints_cvt, $
                              num_constr, $
                              nhard_cvt,  $
                              weights_cvt, $
                              order_cvt,  $
                              knots_cvt,$
                              result1, $
                              result2, $
                              result3, $
                              result4, $
                              result5, $
                              result6, $
                              result7

 result = { DOMAIN_DIM:result1, TARGET_DIM:result2, ORDER:result3, NUM_COEF:result4, NUM_KNOTS:result5, KNOTS:result6, COEF:result7}
   
 RETURN, result
END

