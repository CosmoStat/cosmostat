; $Id: //depot/idl/IDL_70/idl_src/libs/imsl/imsl_6.0/lib/imsl_random_sample.pro#1 $
;
; Copyright (c) 1970-2006, VISUAL NUMERICS Inc. All Rights Reserved.
; This software is confidential information which is proprietary to and a
; trade secret of Visual Numerics Inc.  Use, duplication or disclosure is
; subject to the terms of an appropriate license agreement.
;
FUNCTION imsl_random_sample, nsamp, $                ;INPUT Scalar LONG
                  population, $                      ;INPUT 1-D or 2-D array: floating point
                  Double=double, $                   ;INPUT Scalar ON/OFF flag
                  First_Call=First_Call, $           ;INPUT Scalar ON/OFF flag
                  Additional_Call=Additional_Call, $ ;INPUT Scalar ON/OFF flag
                  Index=Index, $                     ;INPUT/OUTPUT 1-D array: LONG
                  Sample=Sample, $                   ;INPUT/OUTPUT 1-D or 2-D array: floating point
                  npop=npop                          ;INPUT/OUTPUT Scalar LONG

@imsl_init.pro
   ON_ERROR, on_err_action
   ;
   ; Error checking.
   ;   - population must be a scalar, 1-D or 2-D array. 
   ;     Size is [nrow, nvar] or [nvar] ==> nrow=1, or nvar=nrow=1 if a population is scalar.
   ;   - If First_Call is set:
   ;     index and npop are required
   ;     index must be a 1-D array of length nsamp.
   ;   - If Additional_Call is set:
   ;     sample, index and npop are required
   ;     index must be a 1-D array of length nsamp.
   ;     sample must be a 1-D or 2D array of size (nsamp, nvar). 1D if nvar=1.
   ;     
   ;
   nargs = n_params()
   IF (nargs NE 2) THEN $
         message, "Incorrect number of arguments."
   nsamp_cvt = IMSL_LONG(nsamp(0))
   size_pop = IMSL_SIZE(population)
   IF (size_pop(0) GT 2) THEN MESSAGE, "POPULATION must be a scalar, 1-D array, or 2-D array."

   IF (size_pop(0) EQ 0)THEN BEGIN 
	nrow_cvt = IMSL_1
	nvar_cvt = IMSL_1
	pop_dim = 0      
   END
   IF (size_pop(0) EQ 1)THEN BEGIN 
	nrow_cvt = IMSL_1
	nvar_cvt = IMSL_N_ELEMENTS(population)
	pop_dim = 1
   END
   IF (size_pop(0) EQ 2)THEN BEGIN 
	nrow_cvt = size_pop(1)
	nvar_cvt = size_pop(2)
	pop_dim = 2
   END

   IF (KEYWORD_SET(First_Call) EQ TRUE) THEN BEGIN
     IF ((ARG_PRESENT(index) + ARG_PRESENT(npop)) NE 2) THEN $
        MESSAGE, "Keywords INDEX and NPOP are required if FIRST_CALL is set."
     index_cvt = IMSL_LONARR(nsamp_cvt)
     npop_cvt = IMSL_0
     first_set = IMSL_1
   END
   sample_dim = IMSL_0
   IF (KEYWORD_SET(Additional_Call) EQ TRUE) THEN BEGIN
     IF ((KEYWORD_SET(sample) + KEYWORD_SET(index) + KEYWORD_SET(npop)) NE 3) THEN $
        MESSAGE, "Keywords SAMPLE, INDEX and NPOP are required if ADDITIONAL_CALL is set."
     size_index = IMSL_SIZE(index)
     if ((size_index(0) ne 1) OR (N_ELEMENTS(index) ne nsamp_cvt)) THEN $
        MESSAGE, "INDEX must be a 1-D array of length NSAMP."
     size_samp = IMSL_SIZE(sample)
     sample_dim = IMSL_LONG(size_samp(0))
     ; Check SAMPLE
     if ((sample_dim ne 1) AND (sample_dim ne 2)) THEN $
        MESSAGE, "SAMPLE must be a 1-D or 2-D array."
     if ((sample_dim EQ 1)) THEN BEGIN
        if (nvar_cvt NE 1) then MESSAGE, "SAMPLE is not the correct size."
        if (N_ELEMENTS(sample) NE nsamp_cvt) then $
           MESSAGE, "SAMPLE is not the correct size."     
     ENDIF
     if ((sample_dim EQ 2)) THEN BEGIN
        if (nsamp_cvt NE size_samp(1)) then MESSAGE, "SAMPLE is not the correct size."
        if (nvar_cvt NE size_samp(2)) then MESSAGE, "SAMPLE is not the correct size."
     ENDIF
     
     index_cvt = IMSL_LONG(index)
     npop_cvt = IMSL_LONG(npop(0))
     additional_set = IMSL_1
   END
   ;
   ; Decide on what precision to use.
   ; Use the highest precision of the input argument(s).
   ;
   type = TYP_FLOAT
   IF (size_pop(N_ELEMENTS(size_pop)-2) GT TYP_FLOAT) THEN type = TYP_DOUBLE
   IF (KEYWORD_SET(double) EQ true) THEN type = TYP_DOUBLE
   ;
   ; Setup the parameters for the call to the system function
   ;
   ; Floating point arguments and keywords
   IF (type EQ TYP_DOUBLE) THEN BEGIN
      ; 
      ; Input
      if (pop_dim EQ 2) THEN  pop_cvt = DOUBLE(TRANSPOSE(population)) ELSE pop_cvt = DOUBLE(population)  
      if (KEYWORD_SET(SAMPLE)) THEN BEGIN
        if (sample_dim EQ 2) THEN  sample_cvt = DOUBLE(TRANSPOSE(sample)) ELSE sample_cvt = DOUBLE(sample) 
      ENDIF
      result = DBLARR(nvar_cvt, nsamp_cvt)
   END ELSE BEGIN
      ; 
      ; Input
      if (pop_dim EQ 2) THEN  pop_cvt = FLOAT(TRANSPOSE(population)) ELSE pop_cvt = FLOAT(population) 
      if (KEYWORD_SET(SAMPLE)) THEN BEGIN
        if (sample_dim EQ 2) THEN  sample_cvt = FLOAT(TRANSPOSE(sample)) ELSE sample_cvt = FLOAT(sample) 
      ENDIF
      result = FLTARR(nvar_cvt, nsamp_cvt)
   END
   ;
   ; Call the system function.
   ;
    err_status = 0L
  MATHSTAT_305, type,  err_status, pop_cvt, nrow_cvt, nvar_cvt, nsamp_cvt, $
                first_set, $
                additional_set, $
                index_cvt, $
                npop_cvt, $
                sample_cvt, $
                result
   ;
   ; Now copy over all output keywords results.
   ;
   IF (ARG_PRESENT(npop) EQ TRUE) THEN BEGIN
     npop = npop_cvt
     index = index_cvt
     IF (KEYWORD_SET(Additional_Call) EQ TRUE) THEN result = sample_cvt
   END
   if (nsamp_cvt GT 1) THEN result = TRANSPOSE(result)
   ;
   ; Return.
   ;
   RETURN, result
END
   

