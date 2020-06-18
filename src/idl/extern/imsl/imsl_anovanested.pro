; $Id: //depot/idl/IDL_70/idl_src/libs/imsl/imsl_6.0/lib/imsl_anovanested.pro#1 $
;
; Copyright (c) 1970-2006, VISUAL NUMERICS Inc. All Rights Reserved.
; This software is confidential information which is proprietary to and a
; trade secret of Visual Numerics Inc.  Use, duplication or disclosure is
; subject to the terms of an appropriate license agreement.
;
PRO IMSL_L_A3EST, nf, ieq, nl, m, n
    nelem_nl = IMSL_N_ELEMENTS(nl)
    if (ieq EQ 0)  THEN BEGIN
        if (nf GT  2)  THEN BEGIN
            mlast = 1
            m = 1
            l = 0
            n = 1
            for i = 1,  (nf - 1)DO BEGIN
                i_ = i - 1
                for j = 1, n DO BEGIN
                   j_  =  j - 1
                   IF ((l + j_) GE nelem_nl) THEN MESSAGE, "N_LEVELS is not the correct size.", /traceback
                   m = m + nl(l + j_)
                END
                l = l + n
                n = m - mlast
                mlast = m
            END
         END ELSE BEGIN
            n = nl(0)
            m = n + 1
         END
      END ELSE BEGIN
        m = 1
        n  =  1
        IF ((nf-1) GE nelem_nl) THEN MESSAGE, "N_LEVELS is not the correct size.", /traceback
        for i = 1, (nf - 1) DO BEGIN
            i_ = i - 1
            n = n*nl(i_)
            m = m + n
         END
      END
END ; L_A3EST


FUNCTION imsl_anovanested, n_factors, $       ;INPUT Scalar LONG
                     eq_option, $             ;INPUT Scalar LONG
                     n_levels  , $            ;INPUT 1-D array: LONG
                     y, $                     ;INPUT 1-D array: floating point
                     double=double, $         ;INPUT Scalar ON/OFF flag
                     confidence=confidence, $ ;INPUT 1-D Scalar floating point
                     ems=ems, $               ;OUTPUT 1-D array: floating point
                     y_means=y_means, $       ;OUTPUT 1-D array: floating point
                     anova_table=anova_table, $ ;OUTPUT 1-D array: floating point
                     var_comp=var_comp        ;OUTPUT 2-D array: floating point

@imsl_init.pro
   ON_ERROR, on_err_action

   ;
   ; Error checking:
   ; The following checks are performed.
   ;   N_FACTORS must be at least 2.
   ;   EQ_OPTION is treated as a zero/nonzero flag.
   ;   N_LEVELS must be a 1D array, the length depends on the value of EQ_OPTION.
   ;     EQ_OPTION = 0:  set LNL = N_ELEMENTS(N_LEVELS)
   ;     EQ_OPTION = 1:  N_LEVELS is length N_FACTORS.
   ;   Y must be a 1D array of length NOBS.
   ;
   nargs = n_params()
   IF (nargs NE 4) THEN message, "Incorrect number of arguments."

   IF (n_factors(0) LT 2) THEN MESSAGE,  "N_FACTORS mest be at least 2."
   IF (eq_option(0) EQ 0) THEN ieq  =  IMSL_0 ELSE ieq  =  IMSL_1
   size_nl = IMSL_SIZE(n_levels)
   IF (size_nl(0) NE 1) THEN BEGIN
      message, "N_LEVELS must be a 1-D array."
   END
   IF ((ieq EQ 1) AND (N_ELEMENTS(n_levels) NE n_factors)) THEN $
       MESSAGE, "N_LEVELS is not the correct size."
   nl_cvt = IMSL_LONG(n_levels)
   ; Get values of M and N from L_A3EST.  These values are used to help
   ; determine LNL, LNLNF and NOBS.  The logic below is straight out of CNL and FNL.
   IMSL_L_A3EST,  n_factors,  ieq,  nl_cvt,  m,  n
   IF ((m LE 0) OR (n LE 0)) THEN MESSAGE, "An element of N_LEVELS was specified less than one."
   IF (ieq EQ 0) THEN iend  =  m ELSE iend  =  n_factors
   FOR i  =  0,  iend-1 DO BEGIN
      IF (nl_cvt(i) LT 1) THEN MESSAGE, "An element of N_LEVELS was specified less than one."
   END
   LNL  =  m
   lnlnf  =  n
   IF (ieq EQ 0) THEN nobs = TOTAL(nl_cvt((m-n):(m-1))) ELSE nobs = n*nl_cvt(n_factors-1)

   size_y = IMSL_SIZE(y)
   IF (size_y(0) NE 1) THEN BEGIN
      message, "Y must be a 1-D array."
   END
   IF (N_ELEMENTS(y) NE nobs) THEN MESSAGE, "Y is not the correct size."
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
   ; Input LONG keyword(s) and arguments
   nf_cvt  =  IMSL_LONG(n_factors(0))
   ; Output LONG keyword(s)
   ; Cmpute length of Y_MEANS is needed.
   IF (ARG_PRESENT(y_means) EQ true) THEN BEGIN
     IF (ieq NE 0) THEN BEGIN
        nwk = n_levels(0)
        for i=1,n_factors-2 DO  nwk = nwk + nwk*n_levels(i)
        nwk = nwk + 1
     END ELSE BEGIN
        nwk = 1
        for i = 0,lnl-lnlnf-1 DO  nwk = nwk + n_levels(i)
     END
     nymeans = nwk
   END
   ; Floating point arguments and keywords
   IF (type EQ TYP_DOUBLE) THEN BEGIN
      ; Input
      y_cvt = double(y)
      IF (KEYWORD_SET(confidence) EQ TRUE) THEN conf_cvt = DOUBLE(confidence) ELSE conf_cvt = DOUBLE(95.)
      ; Output
      IF (ARG_PRESENT(anova_table) EQ TRUE) THEN at_spc = DBLARR(15)
      IF (ARG_PRESENT(var_comp) EQ TRUE) THEN vc_spc = DBLARR(9, nf_cvt)
      IF (ARG_PRESENT(ems) EQ TRUE) THEN ems_spc = DBLARR((nf_cvt)*((nf_cvt+1)/2))
      IF (ARG_PRESENT(y_means) EQ TRUE) THEN y_means_spc= DBLARR(nymeans)
      result = DOUBLE(0.0)
   END ELSE BEGIN
      ; Input
      y_cvt = FLOAT(y)
      IF (KEYWORD_SET(confidence) EQ TRUE) THEN conf_cvt = FLOAT(confidence) ELSE conf_cvt = FLOAT(95.)
      ; Output
      IF (ARG_PRESENT(anova_table) EQ TRUE) THEN at_spc = FLTARR(15)
      IF (ARG_PRESENT(var_comp) EQ TRUE) THEN vc_spc = FLTARR(9, nf_cvt)
      IF (ARG_PRESENT(ems) EQ TRUE) THEN ems_spc = FLTARR((nf_cvt)*((nf_cvt+1)/2))
      IF (ARG_PRESENT(y_means) EQ TRUE) THEN y_means_spc= FLTARR(nymeans)
      result = float(0.0)
   END
   ;
   ; Call the system function.
   ;
   err_status = 0L
   MATHSTAT_275,type,err_status, nf_cvt, $
                           ieq, $
                           nl_cvt, $
                           y_cvt, $
                           conf_cvt, $
                           at_spc, $
                           vc_spc, $
                           ems_spc, $
                           y_means_spc, $
                           result
   ;
   ; Now copy over all output keywords results.
   ;
   IF (ARG_PRESENT(anova_table) EQ TRUE) THEN anova_table = at_spc
   IF (ARG_PRESENT(var_comp) EQ TRUE) THEN var_comp = TRANSPOSE(vc_spc)
   IF (ARG_PRESENT(ems) EQ TRUE) THEN ems = ems_spc
   IF (ARG_PRESENT(y_means) EQ TRUE) THEN y_means = y_means_spc
   ;
   ; Return.
   ;
   RETURN, result
END









