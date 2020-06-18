; $Id: //depot/idl/IDL_70/idl_src/libs/imsl/imsl_6.0/lib/imsl_csshape.pro#1 $   
;
; Copyright (c) 1970-2006, VISUAL NUMERICS Inc. All Rights Reserved.
; This software is confidential information which is proprietary to and a
; trade secret of Visual Numerics Inc.  Use, duplication or disclosure is
; subject to the terms of an appropriate license agreement.
;
   
FUNCTION imsl_csshape, data1, $               ;INPUT 1-D array: floating point 
                data2, $                   ;INPUT 1-D array: floating point 
                concave=concave, $         ;INPUT Scalar ON/OFF flag
                itmax=itmax, $             ;INPUT Scalar LONG
                double=double              ;INPUT Scalar ON/OFF flag

@imsl_init.pro
   ON_ERROR, on_err_action
   nargs = n_params()
   IF (nargs NE 2) THEN message, 'Incorrect number of arguments.'
   RETURN, imsl_compute_spline(data1, data2, $
                          fcn_idx = IMSL_2, $
                          fcn_name="CSSHAPE", $
                          concave=concave, $
                          itmax=itmax, $  
                          double = double)
END

                   
                   
                   

  
      

  
