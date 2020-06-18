; $Id: //depot/idl/IDL_70/idl_src/libs/imsl/imsl_6.0/lib/imsl_bessy.pro#1 $   
;
; Copyright (c) 1970-2006, VISUAL NUMERICS Inc. All Rights Reserved.
; This software is confidential information which is proprietary to and a
; trade secret of Visual Numerics Inc.  Use, duplication or disclosure is
; subject to the terms of an appropriate license agreement.
;
   
FUNCTION imsl_bessy, order, $                 ;INPUT scalar floating point 
                z, $ ;INPUT scalar or 1-D array: floating point or complex
                double=double, $           ;INPUT Scalar ON/OFF flag
                sequence=sequence          ;INPUT Scalar LONG
@imsl_init.pro
   ON_ERROR, on_err_action

   RETURN, imsl_call_bessel(order, z, $
                       fcn_name = "BESSY", $
                       sequence = sequence, $
                       double = double)
END

