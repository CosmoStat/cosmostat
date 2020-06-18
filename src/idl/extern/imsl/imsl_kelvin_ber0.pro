; $Id: //depot/idl/IDL_70/idl_src/libs/imsl/imsl_6.0/lib/imsl_kelvin_ber0.pro#1 $
;
; Copyright (c) 1970-2006, VISUAL NUMERICS Inc. All Rights Reserved.
; This software is confidential information which is proprietary to and a
; trade secret of Visual Numerics Inc.  Use, duplication or disclosure is
; subject to the terms of an appropriate license agreement.
;
FUNCTION imsl_kelvin_ber0, x, $               ;INPUT 1-D array: floating point
                derivative=derivative, $ ;INPUT Scalar ON/OFF flag
                double=double      ;INPUT Scalar ON/OFF flag
@imsl_init.pro
   ON_ERROR, on_err_action

   RETURN, imsl_call_spcl_fcn1(x, $
                          fcn_name = "KELVIN_BER0", $
                          double = double, $
                          inverse = derivative)
END
