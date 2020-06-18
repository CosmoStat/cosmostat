; $Id: //depot/idl/IDL_70/idl_src/libs/imsl/imsl_6.0/lib/imsl_bsinterp.pro#1 $   
;
; Copyright (c) 1970-2006, VISUAL NUMERICS Inc. All Rights Reserved.
; This software is confidential information which is proprietary to and a
; trade secret of Visual Numerics Inc.  Use, duplication or disclosure is
; subject to the terms of an appropriate license agreement.
;
   
FUNCTION imsl_bsinterp, data1, $              ;INPUT 1-D array: floating point 
                data2, $                   ;INPUT 1-D array: floating point 
                data3, $                   ;INPUT 2-D array: floating point 
                xorder=xorder, $           ;INPUT Scalar LONG
                yorder=yorder, $           ;INPUT Scalar LONG
                double=double, $           ;INPUT Scalar ON/OFF flag
                xknots=xknots, $           ;INPUT 1-D array: floating point
                yknots=yknots              ;INPUT 1-D array: floating point

@imsl_init.pro
   ON_ERROR, on_err_action
   nargs = n_params()
   IF ((nargs NE 2) AND (nargs NE 3)) THEN message, 'Incorrect number of arguments.'
   IF (nargs EQ 2) THEN  $
     RETURN, imsl_compute_spline(data1, data2, /bspline, double = double, $
                            fcn_idx = IMSL_3, $
                            fcn_name = "BSINTERP", $
                            xorder = xorder, yorder = yorder, $
                            xknots = xknots, yknots = yknots) $
   ELSE $
     RETURN, imsl_compute_spline(data1, data2, data3, /bspline, double = double, $
                            fcn_idx = IMSL_4, $
                            fcn_name = "BSINTERP", $
                            xorder = xorder, yorder = yorder, $
                            xknots = xknots, yknots = yknots)
                           
   ; return
END

                   
                   
                   

  
      

  
