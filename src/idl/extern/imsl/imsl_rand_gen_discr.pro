; $Id: //depot/idl/IDL_70/idl_src/libs/imsl/imsl_6.0/lib/imsl_rand_gen_discr.pro#1 $
;
; Copyright (c) 1970-2006, VISUAL NUMERICS Inc. All Rights Reserved.
; This software is confidential information which is proprietary to and a
; trade secret of Visual Numerics Inc.  Use, duplication or disclosure is
; subject to the terms of an appropriate license agreement.
;
FUNCTION imsl_rand_gen_discr, n, $           ;INPUT Scalar LONG
                         imin, $             ;INPUT Scalar LONG
                         nmass, $            ;INPUT Scalar LONG
                         probs, $            ;INPUT 1-D array: floating point
                         table=table, $      ;INPUT Scalar ON/OFF flag
                         double=double       ;INPUT Scalar ON/OFF flag

@imsl_init.pro
   ON_ERROR, on_err_action
   ;
   ; Error checking:
   ;   - N GE 1
   ;   - MASS GT 1
   ;   - PROBS must be a 1-D array
   ;   - If TABLE not set then N-elements(probs) must be EQ (NMASS)
   ;   - If TABLE set then N-elements(probs) must be GE (NMASS+1)
   ;
   nargs = n_params()
   IF (nargs NE 4) THEN $
         message, "Incorrect number of arguments."
   n_cvt = (IMSL_LONG(n))(0) 
   if (n_cvt LT 1) THEN MESSAGE, 'N must be greater than zero.'
   imin_cvt = (IMSL_LONG(imin))(0) 
   nmass_cvt = (IMSL_LONG(nmass))(0) 
   if (nmass_cvt LT 2) THEN MESSAGE, 'NMASS must be greater than one.'

   size_probs = IMSL_LONG(size(probs))
   IF (size_probs(0) NE 1) THEN message, 'PROBS must be a 1-D array.'
   if (KEYWORD_SET(table)) THEN BEGIN
     table_set = IMSL_1
     IF (size_probs(1) LT (nmass_cvt+1)) then MESSAGE, "PROBS is not the correct size."
   END ELSE BEGIN
     IF (size_probs(1) NE (nmass_cvt)) then MESSAGE, "PROBS is not the correct size."
   END
   ;
   ; Decide on what precision to use.
   ; Use the highest precision of the input argument(s).
   ;
   type = TYP_FLOAT
   IF (size_probs(N_ELEMENTS(size_probs)-2) EQ  TYP_DOUBLE) THEN type = TYP_DOUBLE
   IF (KEYWORD_SET(double) EQ true) THEN type = TYP_DOUBLE
   ;
   result = IMSL_LONARR(n_cvt)
   ; Floating point arguments and keywords
   IF (type EQ TYP_DOUBLE) THEN BEGIN
      ; Input
      probs_cvt = DOUBLE(probs)
   END ELSE BEGIN
      ; Input
      probs_cvt = FLOAT(probs)
   END
   ;
   ; Call the system function.
   ;
   err_status = 0L
   MATHSTAT_311, type, err_status, n_cvt, imin_cvt, $
   	nmass_cvt, probs_cvt, table_set, result
   ;
   ; Return.
   ;
   RETURN, result
END
   

