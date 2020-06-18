; $Id: //depot/idl/IDL_70/idl_src/libs/imsl/imsl_6.0/lib/imsl_ranks.pro#1 $      
;
; Copyright (c) 1970-2006, VISUAL NUMERICS Inc. All Rights Reserved.
; This software is confidential information which is proprietary to and a
; trade secret of Visual Numerics Inc.  Use, duplication or disclosure is
; subject to the terms of an appropriate license agreement.
;
FUNCTION imsl_ranks, X, $                    ;INPUT 1-D array: floating point 
           average_tie=average_tie, $     ;INPUT Scalar ON/OFF flag
           blom_scores=blom_scores, $     ;INPUT Scalar ON/OFF flag
           double=double, $               ;INPUT Scalar ON/OFF flag
           exp_norm_scores=exp_norm_scores, $ ;INPUT Scalar ON/OFF flag
           fuzz=fuzz, $                   ;INPUT Scalar floating point
           highest=highest, $             ;INPUT Scalar ON/OFF flag
           lowest=lowest, $               ;INPUT Scalar ON/OFF flag
           random_split=random_split, $   ;INPUT Scalar ON/OFF flag
           ranks=ranks, $                 ;INPUT Scalar ON/OFF flag
           savage_scores=savage_scores, $ ;INPUT Scalar ON/OFF flag
           tukey_scores=tukey_scores, $   ;INPUT Scalar ON/OFF flag
           vdw_scores=vdw_scores          ;INPUT Scalar ON/OFF flag

@imsl_init.pro
   ON_ERROR, on_err_action
   ;
   ; Error checking:
   ;  - X must be either a 1-D array. (set n = length)
   ;  - At most one of AVERAGE_TIE, HIGHEST, LOWEST, or RANDOM_SPLIT
   ;    can be specified.
   ;  - at most one of *_SCORES can be specified.
   ;          
   nargs = n_params()
   IF (nargs NE 1)  THEN message, 'Incorrect number of arguments.'
     
   size_x = IMSL_LONG(size(x))
   IF (size_x(0) NE 1) THEN message, 'X must be a 1-D array.'
   n = IMSL_LONG(size_x(1))
   tmp = KEYWORD_SET(average_tie) + KEYWORD_SET(highest) + $
     KEYWORD_SET(lowest) + KEYWORD_SET(random_split)
   IF (tmp GT 1) THEN message, $
     'At most one of AVERAGE_TIE, HIGHEST, LOWEST, or RANDOM_SPLIT can be set'
   tmp = KEYWORD_SET(blom_scores) + KEYWORD_SET(exp_norm_scores) + $
     KEYWORD_SET(ranks) + KEYWORD_SET(savage_scores) + $
     KEYWORD_SET(tukey_scores) + KEYWORD_SET(vdw_scores)
   IF (tmp GT 1) THEN message, $
     'At most one of RANKS, BLOM_SCORES, TUKEY_SCORES, VDW_SCORES, ' + $
     'EXP_NORM_SCORES, or SAVAGE_SCORES can be set'
   ;
   ; Decide on what precision to use.
   ;
   type = TYP_FLOAT
   IF (size_x(N_ELEMENTS(size_x)-2) EQ  TYP_DOUBLE) THEN type = TYP_DOUBLE
   IF (KEYWORD_SET(double) EQ true) THEN type = TYP_DOUBLE
   ;
   ; Setup the parameters for the call to the system function.
   ;
   ; Floating point arguments and keywords
   ;
   IF (type EQ TYP_DOUBLE) THEN BEGIN
      result = dblarr(n)
      x_cvt = double(x)
      IF (KEYWORD_SET(fuzz)) THEN fuzz_cvt = double(fuzz) $
        ELSE fuzz_cvt = double(0.0)
   END ELSE BEGIN
      result = fltarr(n)
      x_cvt = float(x)
      IF (KEYWORD_SET(fuzz)) THEN fuzz_cvt = float(fuzz) $
        ELSE fuzz_cvt = float(0.0)
   END
   ;
   ; Call the system function.
   ;
   err_status = 0L
   MATHSTAT_183, type, err_status, X_cvt, n, $
                              IMSL_LONG(KEYWORD_SET(average_tie)), $
                              IMSL_LONG(KEYWORD_SET(blom_scores)), $
                              IMSL_LONG(KEYWORD_SET(exp_norm_scores)), $
                              IMSL_LONG(KEYWORD_SET(highest)), $
                              IMSL_LONG(KEYWORD_SET(lowest)), $
                              IMSL_LONG(KEYWORD_SET(random_split)), $
                              IMSL_LONG(KEYWORD_SET(ranks)), $
                              IMSL_LONG(KEYWORD_SET(savage_scores)), $
                              IMSL_LONG(KEYWORD_SET(tukey_scores)), $
                              IMSL_LONG(KEYWORD_SET(vdw_scores)), $
                              fuzz_cvt, $
                              result
   ; Return
   RETURN, result
END
