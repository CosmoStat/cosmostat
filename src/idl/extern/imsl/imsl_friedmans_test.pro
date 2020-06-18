; $Id: //depot/idl/IDL_70/idl_src/libs/imsl/imsl_6.0/lib/imsl_friedmans_test.pro#1 $
;
; Copyright (c) 1970-2006, VISUAL NUMERICS Inc. All Rights Reserved.
; This software is confidential information which is proprietary to and a
; trade secret of Visual Numerics Inc.  Use, duplication or disclosure is
; subject to the terms of an appropriate license agreement.
;
FUNCTION imsl_friedmans_test, y, $            ;INPUT 2-D array: floating point
                     double=double, $         ;INPUT Scalar ON/OFF flag
                     fuzz=fuzz, $             ;INPUT 1-D Scalar floating point
                     alpha=alpha, $           ;INPUT 1-D Scalar floating point
                     diff=diff, $             ;OUTPUT Scalar floating point
                     stats=stats, $             ;OUTPUT 1-D array: floating point
                     sum_rank=sum_rank        ;OUTPUT 1-D array: floating point

@imsl_init.pro
   ON_ERROR, on_err_action

   ;
   ; Error checking:
   ; The following checks are performed.
   ;   Y must be a 2D array   
   ;
   nargs = n_params()
   IF (nargs NE 1) THEN message, "Incorrect number of arguments."
   size_y = IMSL_SIZE(y)
   IF (size_y(0) NE 2) THEN BEGIN
      message, "Y must be a 2-D array."
   END
   nblocks = size_y(1)
   n_treatments = size_y(2)
   ;
   ; Decide on what precision to use.
   ; Use the highest precision of the input argument(s).
   ;
   type = TYP_FLOAT
   IF (size_y(N_ELEMENTS(size_y)-2) GT TYP_FLOAT) THEN type = TYP_DOUBLE
   IF (KEYWORD_SET(double) EQ true) THEN type = TYP_DOUBLE
   ;
   ; Setup the parameters for the call to the system function.
   ;
   ; Input LONG keyword(s)
   ; Output LONG keyword(s)
   ;
   ; Floating point arguments and keywords
   IF (type EQ TYP_DOUBLE) THEN BEGIN
      ; Input
      IF (KEYWORD_SET(fuzz) EQ TRUE) THEN fuzz_cvt = DOUBLE(fuzz) ELSE fuzz_cvt = DOUBLE(0.0)
      IF (KEYWORD_SET(alpha) EQ TRUE) THEN alpha_cvt = DOUBLE(alpha) ELSE alpha_cvt = DOUBLE(0.05)
      ; Output
      IF (ARG_PRESENT(stats) EQ TRUE) THEN stats_spc = DBLARR(6)
      IF (ARG_PRESENT(sum_rank) EQ TRUE) THEN sum_rank_spc = DBLARR(n_treatments)
      IF (ARG_PRESENT(diff) EQ TRUE) THEN diff_spc = DOUBLE(0.0)
      result = DOUBLE(0.0)
   END ELSE BEGIN
      ; Input
      IF (KEYWORD_SET(fuzz) EQ TRUE) THEN fuzz_cvt = FLOAT(fuzz) ELSE fuzz_cvt = FLOAT(0.0)
      IF (KEYWORD_SET(alpha) EQ TRUE) THEN alpha_cvt = FLOAT(alpha) ELSE alpha_cvt = FLOAT(0.05)
      ; Output
      IF (ARG_PRESENT(stats) EQ TRUE) THEN stats_spc = FLTARR(6)
      IF (ARG_PRESENT(sum_rank) EQ TRUE) THEN sum_rank_spc = FLTARR(n_treatments)
      IF (ARG_PRESENT(diff) EQ TRUE) THEN diff_spc = FLOAT(0.0)
      result = FLOAT(0.0)
   END
   y_cvt = imsl_cvt_arr(y, type)
   ;
   ; Call the system function.
   ;
   err_status = 0L
   MATHSTAT_261,type,err_status, y_cvt,  nblocks, n_treatments, $
                           fuzz_cvt, $
                           alpha_cvt, $
                           stats_spc, $
                           sum_rank_spc, $
                           diff_spc, $
                           result
   ;
   ; Now copy over all output keywords results.
   ;
   IF (ARG_PRESENT(diff) EQ TRUE) THEN diff = diff_spc
   IF (ARG_PRESENT(stats) EQ TRUE) THEN stats = stats_spc
   IF (ARG_PRESENT(sum_rank) EQ TRUE) THEN sum_rank = sum_rank_spc
   ;
   ; Return.
   ;
   RETURN, result
END

                   
                   
                   

  
      

  
