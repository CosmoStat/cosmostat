; $Id: //depot/idl/IDL_70/idl_src/libs/imsl/imsl_6.0/lib/imsl_discr_table.pro#1 $
;
; Copyright (c) 1970-2006, VISUAL NUMERICS Inc. All Rights Reserved.
; This software is confidential information which is proprietary to and a
; trade secret of Visual Numerics Inc.  Use, duplication or disclosure is
; subject to the terms of an appropriate license agreement.
;
FUNCTION imsl_discr_table, prf, $              ;INPUT Scalar STRING
                     del, $                    ;INPUT Scalar floating point
                     nndx, $                   ;INPUT Scalar LONG
                     imin, $                   ;INPUT/OUTPUT Scalar LONG
                     nmass, $                  ;INPUT/OUTPUT Scalar LONG
                     cum_probs=cum_probs, $    ;INPUT 1-D array: floating point
                     double=double             ;INPUT Scalar ON/OFF flag
 
@imsl_init.pro
   ON_ERROR, on_err_action
   ;
   ; Error checking:
   ; The following checks are performed.
   ;   - PRF must be a scalar string.
   ;   - If CUM_PROBS is set, then it must be a 1-D array of length nmass.
   ;     - NNDX must be greater than or equal to 1 (If CUM_PROBS set).
   ;
   nargs = n_params()
   IF (nargs LT 5) THEN MESSAGE,  "Incorrect number of arguments."
   size_prf = IMSL_SIZE(prf)
   IF ((N_ELEMENTS(prf) NE 1) OR (size_prf(N_ELEMENTS(size_prf)-2) NE 7)) THEN $
     message, 'PRF must be a scalar string.'
   ;
   nndx_cvt = IMSL_LONG(nndx(0))
   imin_cvt = IMSL_LONG(imin(0))
   nmass_cvt = IMSL_LONG(nmass(0))   
   ;
   IF (KEYWORD_SET(cum_probs)) THEN BEGIN
      size_cum_probs  = IMSL_SIZE(cum_probs)
      IF (size_cum_probs(0) NE 1) THEN BEGIN
         message, "CUM_PROBS must be a 1-D array."
      END
      IF (N_ELEMENTS(cum_probs) NE nmass_cvt) THEN MESSAGE, "CUM_PROBS is not the correct length."
      IF (nndx_cvt LT 0) THEN MESSAGE, "NNDX must be greater than or equal to 1."
   END  
   ;   
   size_del = IMSL_LONG(size(del))
   ;
   ; Decide on what precision to use.
   ; Use the highest precision of the input argument(s).
   ;
   type = TYP_FLOAT
   IF (size_del(N_ELEMENTS(size_del)-2) EQ  TYP_DOUBLE) THEN type = TYP_DOUBLE
   IF (KEYWORD_SET(double) EQ true) THEN type = TYP_DOUBLE
   ;
   if (KEYWORD_SET(cum_probs)) THEN index_cvt = IMSL_1 
   ; Floating point arguments and keywords
   IF (type EQ TYP_DOUBLE) THEN BEGIN
      ; Input
      del_cvt = DOUBLE(del(0))
      if (KEYWORD_SET(cum_probs)) THEN BEGIN
        result = DBLARR(N_ELEMENTS(cum_probs)+nndx_cvt)
	result(0:N_ELEMENTS(cum_probs)-1) = DOUBLE(cum_probs)
      END ELSE BEGIN
        result = DOUBLE(0)
      END
   END ELSE BEGIN
      ; Input
      del_cvt = FLOAT(del(0))
      if (KEYWORD_SET(cum_probs)) THEN BEGIN
        result = FLTARR(N_ELEMENTS(cum_probs)+nndx_cvt)
	result(0:N_ELEMENTS(cum_probs)-1) = FLOAT(cum_probs)
      END ELSE BEGIN
        result = FLOAT(0)
      END
   END
   ;
   ; Call the system function.
   ;
   err_status = 0L
   MATHSTAT_309,  type,  err_status,  $
                                      prf,   $
                                      del_cvt,   $
                                      nndx_cvt,   $
                                      imin_cvt,   $
                                      nmass_cvt,   $
                                      index_cvt,   $
                                      result, $
                                      IMSL_N_ELEMENTS(result)


   imin = imin_cvt
   nmass = nmass_cvt
   ;
   ; Return.
   ;
   RETURN, result
END

                   
                   
                   

  
      

  
