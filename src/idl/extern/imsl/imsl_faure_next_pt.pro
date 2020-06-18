; $Id: //depot/idl/IDL_70/idl_src/libs/imsl/imsl_6.0/lib/imsl_faure_next_pt.pro#1 $
;
; Copyright (c) 1970-2006, VISUAL NUMERICS Inc. All Rights Reserved.
; This software is confidential information which is proprietary to and a
; trade secret of Visual Numerics Inc.  Use, duplication or disclosure is
; subject to the terms of an appropriate license agreement.
;
FUNCTION imsl_faure_next_pt, npts, $                  ;INPUT Scalar LONG
                  state, $                            ;INPUT Scalar LONG
                  skip=skip, $                        ;OUTPUT Scalar LONG
                  Double=DOUBLE                       ;INPUT Scalar ON/OFF flag
@imsl_init.pro
   ON_ERROR, on_err_action
   ;
   ; Error checking.
   ; NPTS > 0
   ; Do some sanity checks on the state structure.
   ;
   nargs = n_params()
   IF (nargs NE 2) THEN $
     MESSAGE, "Incorrect number of arguments."
   npts_cvt = IMSL_LONG(npts(0))
   IF (npts_cvt LT 1) THEN MESSAGE, "NPTS must be positive."
   ; Always get Skip, only return it if needed
   skip_spc = IMSL_0
   ;
   ; Decide on what precision to use.
   ;
   type = TYP_FLOAT
   IF (KEYWORD_SET(DOUBLE) EQ true) THEN type = TYP_DOUBLE
   ndim = state.dim
   IF (type EQ TYP_DOUBLE) THEN result = DBLARR(ndim, npts_cvt) ELSE result = FLTARR(ndim, npts_cvt)
   ;
   ; Call the system function.
   ;
   nskip_cvt = state.NSKIP
   y_cvt = state.Y
   dim_cvt = state.DIM
   maxdigits_cvt = state.MAXDIGITS
   base_cvt = state.BASE
   digitsn_cvt = state.DIGITSN
   c_cvt = state.C
   power_cvt = state.POWER
   scale_cvt = state.SCALE
   err_status = 0L
   MATHSTAT_315, type, err_status, $
                npts_cvt, $
                skip_spc, $
                NSKIP_cvt, $
                Y_cvt, $
                DIM_cvt, $
                MAXDIGITS_cvt, $
                BASE_cvt, $
                DIGITSN_cvt, $
                C_cvt, $
                POWER_cvt, $
                SCALE_cvt, $
                result
   ;
   ; Return.
   ;
   IF (ARG_PRESENT(skip) EQ TRUE) THEN skip = skip_spc
   ;
   ; Update the stat structure for future calls.
   state = {NSKIP:nskip_cvt, $
               Y:y_cvt, DIM:dim_cvt, $
               MAXDIGITS:maxDigits_cvt, $
               BASE:base_cvt, $
               DIGITSN:digitsN_cvt, $
               C:c_cvt, $
               POWER:power_cvt, $
               SCALE:scale_cvt}
   
   RETURN, TRANSPOSE(result)
END
   

