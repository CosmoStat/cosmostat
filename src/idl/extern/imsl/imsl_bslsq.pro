; $Id: //depot/idl/IDL_70/idl_src/libs/imsl/imsl_6.0/lib/imsl_bslsq.pro#1 $   
;
; Copyright (c) 1970-2006, VISUAL NUMERICS Inc. All Rights Reserved.
; This software is confidential information which is proprietary to and a
; trade secret of Visual Numerics Inc.  Use, duplication or disclosure is
; subject to the terms of an appropriate license agreement.
;
   
FUNCTION imsl_bslsq, data1, $              ;INPUT 1-D array: floating point 
                data2, $                   ;INPUT 1-D array: floating point 
                data3, $                   ;INPUT 2-D array: floating point 
                data4, $
                data5, $
                sse=sse, $                 ;INPUT Scalar floating point
                optimize=optimize, $       ;INPUT Scalar ON/OFF flag
                xweights=xweights, $       ;INPUT 1-D array: floating point
                yweights=yweights, $       ;INPUT 1-D array: floating point
                xorder=xorder, $           ;INPUT Scalar LONG
                yorder=yorder, $           ;INPUT Scalar LONG
                double=double, $           ;INPUT Scalar ON/OFF flag
                xknots=xknots, $           ;INPUT 1-D array: floating point
                yknots=yknots              ;INPUT 1-D array: floating point

@imsl_init.pro
   ON_ERROR, on_err_action
   nargs = n_params()
   IF ((nargs NE 5) AND (nargs NE 3)) THEN message, 'Incorrect number of arguments.'
   IF (nargs EQ 3) THEN  $
     RETURN, imsl_compute_spline(data1, data2, $
                            fcn_idx = IMSL_LONG(5), $
                            fcn_name = "BSLSQ", $
                            sse = sse, $
                            optimize=optimize, $
                            xweights=xweights, $
                            xspace_dim = IMSL_LONG(data3), $
                            /bspline, double = double, $
                            xorder = xorder, yorder = yorder, $
                            xknots = xknots, yknots = yknots) $
   ELSE $
     RETURN, imsl_compute_spline(data1, data2, data3, $
                            fcn_idx = IMSL_LONG(6), $
                            fcn_name = "BSLSQ", $
                            sse = sse, $
                            optimize=optimize, $
                            xweights=xweights, $
                            yweights=yweights, $
                            xspace_dim = IMSL_LONG(data4), $
                            yspace_dim = IMSL_LONG(data5), $
                            /bspline, double = double, $
                            xorder = xorder, yorder = yorder, $
                            xknots = xknots, yknots = yknots)
                           
   ; return
END

                   
                   
                   

  
      

  
