; $Id: //depot/idl/IDL_70/idl_src/libs/imsl/imsl_6.0/lib/imsl_bessi_exp.pro#1 $
;
; Copyright (c) 1970-2006, VISUAL NUMERICS Inc. All Rights Reserved.
; This software is confidential information which is proprietary to and a
; trade secret of Visual Numerics Inc.  Use, duplication or disclosure is
; subject to the terms of an appropriate license agreement.
;
FUNCTION imsl_bessi_exp, order, $       ;INPUT scalar int, either 0 or 1. 
		x, $               ;INPUT 1-D array: floating point
                double=double      ;INPUT Scalar ON/OFF flag
@imsl_init.pro
   ON_ERROR, on_err_action

   order_cvt = IMSL_LONG(order(0))

   if ((order_cvt NE 0) and (order_cvt NE 1)) then $
      message, "ORDER must be either zero or one."

   if (order_cvt EQ 0) then begin
      result =  imsl_call_spcl_fcn1(x, $
                          fcn_name = "BESSI0_EXP", $
                          double = double)
   endif
   if (order_cvt EQ 1) then begin
      result =  imsl_call_spcl_fcn1(x, $
                          fcn_name = "BESSI1_EXP", $
                          double = double)
   endif
   return, result
END
