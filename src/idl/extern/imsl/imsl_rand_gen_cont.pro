; $Id: //depot/idl/IDL_70/idl_src/libs/imsl/imsl_6.0/lib/imsl_rand_gen_cont.pro#1 $
;
; Copyright (c) 1970-2006, VISUAL NUMERICS Inc. All Rights Reserved.
; This software is confidential information which is proprietary to and a
; trade secret of Visual Numerics Inc.  Use, duplication or disclosure is
; subject to the terms of an appropriate license agreement.
;
FUNCTION imsl_rand_gen_cont, n, $          ;INPUT Scalar LONG
                     table, $              ;INPUT 2-D array: floating point
                   double=double            ;INPUT Scalar ON/OFF flag

@imsl_init.pro
   ON_ERROR, on_err_action
   ;
   ; Error checking:
   ;   - N GE 1
   ;   - TABLE must be NDATA by 5.
   ;
   nargs = n_params()
   IF (nargs NE 2) THEN $
         message, "Incorrect number of arguments."
   n_cvt = (IMSL_LONG(n))(0) 
   if (n_cvt LT 1) THEN MESSAGE, 'N must be greater than zero.'

   size_table = IMSL_LONG(size(table))
   IF (size_table(0) NE 2) THEN message, 'TABLE must be a 2-D array.'
   ndata_cvt = size_table(1) 
   if (size_table(2) NE 5) then MESSAGE, "TABLE is not the correct size."
   ;
   ; Decide on what precision to use.
   ; Use the highest precision of the input argument(s).
   ;
   type = TYP_FLOAT
   IF (size_table(N_ELEMENTS(size_table)-2) EQ  TYP_DOUBLE) THEN type = TYP_DOUBLE
   IF (KEYWORD_SET(double) EQ true) THEN type = TYP_DOUBLE
   ;
   ; Floating point arguments and keywords
   IF (type EQ TYP_DOUBLE) THEN BEGIN
      ; Input
      table_cvt = DOUBLE(TRANSPOSE(table))
      result = dblarr(n_cvt)
   END ELSE BEGIN
      ; Input
      table_cvt = FLOAT(TRANSPOSE(table))
      result = fltarr(n_cvt)
   END
   ;
   ; Call the system function.
   ;
   err_status = 0L
   MATHSTAT_308, type, err_status, n_cvt, ndata_cvt, table_cvt, result
   ;
   ; Return.
   ;
   RETURN, result
END
   

