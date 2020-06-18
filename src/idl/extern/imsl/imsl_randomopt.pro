; $Id: //depot/idl/IDL_70/idl_src/libs/imsl/imsl_6.0/lib/imsl_randomopt.pro#1 $   
;
; Copyright (c) 1970-2006, VISUAL NUMERICS Inc. All Rights Reserved.
; This software is confidential information which is proprietary to and a
; trade secret of Visual Numerics Inc.  Use, duplication or disclosure is
; subject to the terms of an appropriate license agreement.
;
PRO imsl_randomopt, gen_option=gen_option, $              ;INPUT Scalar LONG
                    get_seed=get_seed, $             ;INPUT Scalar LONG
                    set_seed=set_seed, $             ;INPUT Scalar LONG
                    substream_seed=substream_seed, $ ;OUTPUT Scalar LONG
                    current_option=current_option    ;INPUT Scalar LONG
@imsl_init.pro
   ON_ERROR, on_err_action
   ;
   ; Error checking: 
   ;  - If SUBSTREAM_SEED is set, then GET_SEED must also be set.
   ;
   if (KEYWORD_SET(substream_seed) EQ TRUE) THEN BEGIN
      if (NOT ARG_PRESENT(GET_SEED)) THEN MESSAGE, $
         "If SUBSTREAM_SEED is set, then GET_SEED must also be used."
   END

   ; Setup the parameters for the call to the system function.
   ;
   ; Input LONG keyword(s)
   IF (KEYWORD_SET(gen_option)) THEN gen_option_cvt = IMSL_LONG(gen_option(0))
   IF (KEYWORD_SET(set_seed)) THEN set_seed_cvt = IMSL_LONG(set_seed(0))
   IF (KEYWORD_SET(substream_seed)) THEN sub_seed_cvt = IMSL_LONG(substream_seed(0))
   ;
   ; Output LONG keyword(s)
   IF (arg_present(get_seed)) THEN get_seed_spc = IMSL_0;
   IF (arg_present(current_option)) THEN current_opt_spc = IMSL_0;
   ;
   ; Call the system function.
   ;
   err_status = 0L
   type = TYP_MEMINT
   MATHSTAT_182, type, err_status, gen_option_cvt, set_seed_cvt, get_seed_spc, current_opt_spc, $
                 sub_seed_cvt

   IF (ARG_PRESENT(get_seed)) THEN get_seed=get_seed_spc
   IF (ARG_PRESENT(current_option)) THEN current_option=current_opt_spc
   ; return
   RETURN
END

                   
                   
                   

  
      

  
