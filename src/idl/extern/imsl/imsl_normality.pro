; $Id: //depot/idl/IDL_70/idl_src/libs/imsl/imsl_6.0/lib/imsl_normality.pro#1 $
;
; Copyright (c) 1970-2006, VISUAL NUMERICS Inc. All Rights Reserved.
; This software is confidential information which is proprietary to and a
; trade secret of Visual Numerics Inc.  Use, duplication or disclosure is
; subject to the terms of an appropriate license agreement.
;
FUNCTION imsl_normality, x, $                      ;INPUT 1-D array: floating point
                    double=double, $            ;INPUT Scalar ON/OFF flag
                    ncat=ncat, $                ;INPUT Scalar LONG
                    df=df, $                    ;OUTPUT Scalar floating point
                    chisq=chisq, $              ;OUTPUT Scalar floating point
                    Lilliefors=lilliefors, $    ;OUTPUT Scalar floating point
                    Shapiro_wilk=shapiro_wilk   ;OUTPUT Scalar floating point

@imsl_init.pro
   ON_ERROR, on_err_action
   ;
   ; Error checking.
   ; - Make sure the input argument is a simple array.
   ; - NCAT, DF, and CHISQ must be used together.
   ;
   nargs = n_params()
   IF (nargs EQ 1) THEN BEGIN
      size_x = IMSL_SIZE(x)
      IF (size_x(0) NE 1) THEN BEGIN
         message, "X must be a 1-D array."
      END
   END ELSE MESSAGE,  "Incorrect number of arguments."
   chsq_nargs = ARG_PRESENT(ncat) + ARG_PRESENT(df) + ARG_PRESENT(chisq)
   IF ((chsq_nargs EQ 1) OR (chsq_nargs EQ 2)) THEN $
     MESSAGE,  "Keywords NCAT, DF and CHISQ must be used together."
   IF (ARG_PRESENT(ncat)) THEN BEGIN 
      IF (N_ELEMENTS(ncat) EQ 0) THEN MESSAGE,  "NCAT must be defined."
   END
   ;
   ; Decide on what precision to use.
   ; Use the highest precision of the input argument(s).
   ;
   type = TYP_FLOAT
   IF (size_x(N_ELEMENTS(size_x)-2) GT TYP_FLOAT) THEN type = TYP_DOUBLE
   IF (KEYWORD_SET(double) EQ true) THEN type = TYP_DOUBLE
   ;
   ; Setup the parameters for the call to the system function.
   ;
   ; Input LONG argument(s)
   IF (ARG_PRESENT(ncat)) THEN ncat_cvt = IMSL_LONG(ncat(0))
   nobs = size_x(N_ELEMENTS(size_x)-1)
   IF (type EQ TYP_DOUBLE) THEN BEGIN
      result = double(0.0)
      x_cvt = double(x)
      ;
      ; Output keywords.
      ;
      IF (ARG_PRESENT(Lilliefors) EQ TRUE) THEN $
        lilliefors_spc = double(0.0)
      IF (ARG_PRESENT(Shapiro_wilk) EQ TRUE) THEN $
        shapiro_wlk_spc = double(0.0)
      IF (ARG_PRESENT(df) EQ TRUE) THEN $
        df_spc = double(0.0)
      IF (ARG_PRESENT(chisq) EQ TRUE) THEN $
        chisq_spc = double(0.0)
   END ELSE BEGIN
      result = float(0.0)
      x_cvt = float(x)
      ;
      ; Output keywords.
      ;
      IF (ARG_PRESENT(Lilliefors) EQ TRUE) THEN $
        lilliefors_spc = float(0.0)
      IF (ARG_PRESENT(Shapiro_wilk) EQ TRUE) THEN $
        shapiro_wlk_spc = float(0.0)
      IF (ARG_PRESENT(df) EQ TRUE) THEN $
        df_spc = float(0.0)
      IF (ARG_PRESENT(chisq) EQ TRUE) THEN $
        chisq_spc = float(0.0)
   END
   ;
   ; Call the system function.
   ;
   err_status = 0L
   MATHSTAT_170, type, err_status, x_cvt, nobs, $
                           lilliefors_spc, shapiro_wlk_spc, $
                           ncat_cvt,  df_spc,  chisq_spc,  result
   ;
   ; Now copy over all output keywords results.
   ;
   IF (ARG_PRESENT(Lilliefors) EQ TRUE) THEN $
     lilliefors=lilliefors_spc
   IF (ARG_PRESENT(Shapiro_wilk) EQ TRUE) THEN $
     shapiro_wilk=shapiro_wlk_spc
   IF (ARG_PRESENT(df) EQ TRUE) THEN $
     df=df_spc
   IF (ARG_PRESENT(chisq) EQ TRUE) THEN $
     chisq=chisq_spc

   ;
   ; Return.
   ;
   RETURN, result
END

                   
                   
                   

  
      

  
