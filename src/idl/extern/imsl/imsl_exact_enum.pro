; $Id: //depot/idl/IDL_70/idl_src/libs/imsl/imsl_6.0/lib/imsl_exact_enum.pro#1 $
;
; Copyright (c) 1970-2006, VISUAL NUMERICS Inc. All Rights Reserved.
; This software is confidential information which is proprietary to and a
; trade secret of Visual Numerics Inc.  Use, duplication or disclosure is
; subject to the terms of an appropriate license agreement.
;
FUNCTION imsl_exact_enum, table, $        ;INPUT 2-D array: floating point 
           p_value=p_value,  $            ;OUTPUT Scalar floating point
           prob_table=prob_table,  $      ;OUTPUT Scalar floating point
           error_chk=error_chk,  $        ;OUTPUT Scalar floating point
           double=double                  ;INPUT Scalar ON/OFF flag
           
@imsl_init.pro
   ON_ERROR, on_err_action
   ;
   ; Error checking:
   ;  - TABLE must be either a 1-D or 2-D array. (m x n)
   ;          
   nargs = n_params()
   IF (nargs NE 1)  THEN message, 'Incorrect number of arguments.'
     
   size_table = IMSL_LONG(size(table))
   IF ((size_table(0) NE 1) AND (size_table(0) NE 2)) THEN $
     message, 'TABLE must be a 1-D or 2-D array.'
   m = IMSL_LONG(size_table(1))
   IF (size_table(0) EQ 2) THEN n = IMSL_LONG(size_table(2)) ELSE n = IMSL_1
   ;
   ; Decide on what precision to use.
   ;
   type = TYP_FLOAT
   IF (size_table(N_ELEMENTS(size_table)-2) EQ  TYP_DOUBLE) THEN type = TYP_DOUBLE
   IF (KEYWORD_SET(double) EQ true) THEN type = TYP_DOUBLE
   ;
   ; Setup the parameters for the call to the system function.
   ;
   ; Floating point arguments and keywords
   ;
   IF (type EQ TYP_DOUBLE) THEN BEGIN
      result = DOUBLE(0.)
      p_value_spc = double(0.)
      prob_table_spc = double(0.)
      error_chk_spc = double(0.)
   END ELSE BEGIN
      result = FLOAT(0.)
      p_value_spc = float(0.)
      prob_table_spc = float(0.)
      error_chk_spc = float(0.)
   END
   table_cvt = imsl_cvt_arr(table, type)
   ;
   ; Call the system function.
   ;
   err_status = 0L
   MATHSTAT_253,    type,    err_status,    table_cvt,    m,    n,  $
     p_value_spc, prob_table_spc,  error_chk_spc, $
     result
   ;
   p_value = p_value_spc
   prob_table = prob_table_spc
   error_chk = error_chk_spc
   ;
   ; Return
   RETURN, result
END
