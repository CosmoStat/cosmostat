; $Id: //depot/idl/IDL_70/idl_src/libs/imsl/imsl_6.0/lib/imsl_randomness_test.pro#1 $
;
; Copyright (c) 1970-2006, VISUAL NUMERICS Inc. All Rights Reserved.
; This software is confidential information which is proprietary to and a
; trade secret of Visual Numerics Inc.  Use, duplication or disclosure is
; subject to the terms of an appropriate license agreement.
;
FUNCTION imsl_randomness_test, x, $           ;INPUT 1-D array: floating point
                     n_run, $                 ;INPUT Scalar LONG
                     double=double, $         ;INPUT Scalar ON/OFF flag
                     runs_counts=runs_counts, $ ;OUTPUT 1-D Scalar floating point
                     covariances=covariances, $       ;OUTPUT 2-D Scalar floating point
                     pairs_counts=pairs_counts, $     ;OUTPUT 2-D Scalar floating point
                     dsquare_counts=dsquare_counts, $ ;OUTPUT 1-D Scalar floating point
                     dcube_counts=dcube_counts, $     ;OUTPUT 2-D Scalar floating point
                     chisq=chisq, $                   ;OUTPUT Scalar floating point
                     df=df, $                         ;OUTPUT Scalar floating point
                     expect=expect, $                 ;OUTPUT Scalar floating point
                     runs_expect=runs_expect, $       ;OUTPUT 1-D Scalar floating point
                     pairs_lag=pairs_lag              ;INPUT Scalar LONG

@imsl_init.pro
   ON_ERROR, on_err_action

   ;
   ; Error checking:
   ; The following checks are performed.
   ;   X must be a 1D array
   ;   RUNS_COUNTS and COVARIANCES must be used together.
   ;   PAIRS_LAG and PAIRS_COUNTS must be used together.
   ;   RUNS_EXPECT is valid only if RUNS_COUNTS is also used.
   ;   EXPECT and RUNS_COUNTS are mutually exclusive.
   ;   only one method can be specified.
   ; 
   ;
   nargs = n_params()
   IF (nargs NE 2) THEN message, "Incorrect number of arguments."
   size_x = IMSL_SIZE(x)
   IF (size_x(0) NE 1) THEN BEGIN
      message, "X must be a 1-D array."
   END
   nobs  =  IMSL_N_ELEMENTS(x)
   IF ((ARG_PRESENT(runs_counts) + ARG_PRESENT(covariances)) EQ 1) then $
     MESSAGE, "RUNS_COUNTS and COVARIANCES must be used together."
   IF ((ARG_PRESENT(PAIRS_LAG) + ARG_PRESENT(PAIRS_COUNTS)) EQ 1) then $
     MESSAGE, "PAIRS_LAG and PAIRS_COUNTS must be used together."
   n_methods = (ARG_PRESENT(runs_counts) + $
      ARG_PRESENT(PAIRS_LAG) + $
      ARG_PRESENT(dsquare_counts) + $
      ARG_PRESENT(dcube_counts))
   IF (n_methods GT 1) THEN $
      MESSAGE, "The keyword supplied imply conflicting test to be performed."
   IF (n_methods EQ 0) THEN method   =   IMSL_1
   IF (ARG_PRESENT(runs_counts)) THEN method = IMSL_1
   IF (ARG_PRESENT(pairs_counts)) THEN method = IMSL_2
   IF (ARG_PRESENT(dsquare_counts)) THEN method = IMSL_3
   IF (ARG_PRESENT(dcube_counts)) THEN method = IMSL_4
   IF (ARG_PRESENT(runs_expect) AND (method NE 1)) THEN MESSAGE, $
      "RUNS_EXPECT is not valid unless a runs test is being performed."
   IF (ARG_PRESENT(expect) AND (method EQ 1)) THEN MESSAGE, $
      "EXPECT is not valid if a runs test is being performed."
   
   ; Decide on what precision to use.
   ; Use the highest precision of the input argument(s).
   ;
   type = TYP_FLOAT
   IF (size_x(N_ELEMENTS(size_x)-2) GT TYP_FLOAT) THEN type = TYP_DOUBLE
   IF (KEYWORD_SET(double) EQ true) THEN type = TYP_DOUBLE
   ;
   ; Setup the parameters for the call to the system function.
   ;
   ; Input LONG keyword(s)
   IF (KEYWORD_SET(pairs_lag)) THEN pairs_lag_cvt  =  IMSL_LONG(pairs_lag(0))
   n_run_cvt = IMSL_LONG(n_run(0))
   ; Output LONG keyword(s)
   ;
   ; Floating point arguments and keywords
   IF (type EQ TYP_DOUBLE) THEN BEGIN
      ; Input
      x_cvt  =  DOUBLE(x)
      ; Output
      IF (ARG_PRESENT(runs_counts) EQ TRUE) THEN rc_spc = DBLARR(n_run_cvt)
      IF (ARG_PRESENT(covariances) EQ TRUE) THEN cov_spc = DBLARR(n_run_cvt, n_run_cvt)
      IF (ARG_PRESENT(pairs_counts) EQ TRUE) THEN pc_spc = DBLARR(n_run_cvt, n_run_cvt)
      IF (ARG_PRESENT(dsquare_counts) EQ TRUE) THEN dsc_spc = DBLARR(n_run_cvt)
      IF (ARG_PRESENT(dcube_counts) EQ TRUE) THEN dcc_spc = DBLARR(n_run_cvt, n_run_cvt, n_run_cvt)
      IF (ARG_PRESENT(runs_expect) EQ TRUE) THEN re_spc = DBLARR(n_run_cvt)
      IF (ARG_PRESENT(chisq) EQ TRUE) THEN chisq_spc  =  DOUBLE(0)
      IF (ARG_PRESENT(df) EQ TRUE) THEN df_spc  =  DOUBLE(0)
      IF (ARG_PRESENT(expect) EQ TRUE) THEN expect_spc  =  DOUBLE(0)
      result = DOUBLE(0)
   END ELSE BEGIN
      ; Input
      x_cvt  =  FLOAT(x)
      ; Output
      IF (ARG_PRESENT(runs_counts) EQ TRUE) THEN rc_spc = FLTARR(n_run_cvt)
      IF (ARG_PRESENT(covariances) EQ TRUE) THEN cov_spc = FLTARR(n_run_cvt, n_run_cvt)
      IF (ARG_PRESENT(pairs_counts) EQ TRUE) THEN pc_spc = FLTARR(n_run_cvt, n_run_cvt)
      IF (ARG_PRESENT(dsquare_counts) EQ TRUE) THEN dsc_spc = FLTARR(n_run_cvt)
      IF (ARG_PRESENT(dcube_counts) EQ TRUE) THEN dcc_spc = FLTARR(n_run_cvt, n_run_cvt, n_run_cvt)
      IF (ARG_PRESENT(runs_expect) EQ TRUE) THEN re_spc = FLTARR(n_run_cvt)
      IF (ARG_PRESENT(chisq) EQ TRUE) THEN chisq_spc  =  FLOAT(0)
      IF (ARG_PRESENT(df) EQ TRUE) THEN df_spc  =  FLOAT(0)
      IF (ARG_PRESENT(expect) EQ TRUE) THEN expect_spc  =  FLOAT(0)
      result = FLOAT(0)
   END
   ;
   ; Call the system function.
   ;
   err_status = 0L
   MATHSTAT_277,type,err_status, x_cvt,  nobs, $
                           n_run_cvt, $
                           method, $
                           rc_spc, $
                           cov_spc, $
                           pairs_lag_cvt, $
                           pc_spc, $
                           dsc_spc, $
                           dcc_spc, $
                           re_spc, $
                           chisq_spc, $
                           df_spc, $
                           expect_spc, $
                           result
   ;
   ; Now copy over all output keywords results.
   ;
   IF (ARG_PRESENT(runs_counts) EQ TRUE) THEN runs_counts = rc_spc
   IF (ARG_PRESENT(covariances) EQ TRUE) THEN covariances = TRANSPOSE(cov_spc)
   IF (ARG_PRESENT(pairs_counts) EQ TRUE) THEN pairs_counts = TRANSPOSE(pc_spc)
   IF (ARG_PRESENT(dsquare_counts) EQ TRUE) THEN dsquare_counts = dsc_spc
   IF (ARG_PRESENT(dcube_counts) EQ TRUE) THEN dcube_counts = TRANSPOSE(dcc_spc)
   IF (ARG_PRESENT(runs_expect) EQ TRUE) THEN runs_expect = re_spc
   IF (ARG_PRESENT(chisq) EQ TRUE) THEN chisq = chisq_spc 
   IF (ARG_PRESENT(df) EQ TRUE) THEN df = df_spc 
   IF (ARG_PRESENT(expect) EQ TRUE) THEN expect = expect_spc 
   IF (ARG_PRESENT(dcube_counts) EQ TRUE) THEN BEGIN
      IF (type EQ TYP_DOUBLE) THEN dcube_counts = DBLARR(n_run_cvt,n_run_cvt,n_run_cvt) $
        ELSE  dcube_counts = FLTARR(n_run_cvt,n_run_cvt,n_run_cvt) 
      FOR i = 0, n_run_cvt-1 DO FOR j   =   0,   n_run_cvt-1 DO $
         dcube_counts(*,  j,  i) = dcc_spc(i, j, *)
   END
   ;
   ; Return.
   ;
   RETURN, result
END

                   
                   
                   

  
      

  
