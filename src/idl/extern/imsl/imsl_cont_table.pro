; $Id: //depot/idl/IDL_70/idl_src/libs/imsl/imsl_6.0/lib/imsl_cont_table.pro#1 $
;
; Copyright (c) 1970-2006, VISUAL NUMERICS Inc. All Rights Reserved.
; This software is confidential information which is proprietary to and a
; trade secret of Visual Numerics Inc.  Use, duplication or disclosure is
; subject to the terms of an appropriate license agreement.
;
PRO imsl_cont_table, f, $                 ;INPUT Scalar STRING
                     iopt, $              ;INPUT Scalar LONG
                     ndata, $             ;INPUT Scalar LONG
                     table, $             ;INPUT/OUTPUT 2-D array: floating point
                     double=double        ;INPUT Scalar ON/OFF flag
 
@imsl_init.pro
   ON_ERROR, on_err_action
   ;
   ; Error checking:
   ; The following checks are performed.
   ;   - F must be a scalar string.
   ;   - IOPT either 0 or 1.
   ;   - NDATA must be GE 4.
   ;   - TABLE must be NDATA by 5.
   ;
   nargs = n_params()
   IF (nargs LT 4) THEN MESSAGE,  "Incorrect number of arguments."
   size_f = IMSL_SIZE(f)
   IF ((N_ELEMENTS(f) NE 1) OR (size_f(N_ELEMENTS(size_f)-2) NE 7)) THEN $
     message, 'F must be a scalar string.'
   iopt_cvt = IMSL_LONG(iopt(0))
   IF ((iopt_cvt NE 0) AND (iopt_cvt NE 1)) THEN MESSAGE,  "IOPT must be either 0 or 1."
   ndata_cvt = IMSL_LONG(ndata(0))
   IF (ndata_cvt LT 4) THEN MESSAGE,  "NDATA must be atleast 4."
   size_table = IMSL_LONG(size(table))
   IF (size_table(0) NE 2) THEN message, 'TABLE must be a 2-D array.'
   if ((size_table(1) NE ndata_cvt) OR (size_table(2) NE 5)) then MESSAGE, "TABLE is not the correct size."
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
   END ELSE BEGIN
      ; Input
      table_cvt = FLOAT(TRANSPOSE(table))
   END
   ;
   ; Call the system function.
   ;
   ne_spc = IMSL_0
   err_status = 0L
   MATHSTAT_307,  type,  err_status,  $
                                      f,   $
                                      iopt_cvt,   $
                                      ndata_cvt,   $
                                      table_cvt

   table = transpose(table_cvt)
   ;
   ; Return.
   ;
   RETURN
END

                   
                   
                   

  
      

  
