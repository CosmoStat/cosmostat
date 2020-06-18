; $Id: //depot/idl/IDL_70/idl_src/libs/imsl/imsl_6.0/lib/imsl_csinterp.pro#1 $   
;
; Copyright (c) 1970-2006, VISUAL NUMERICS Inc. All Rights Reserved.
; This software is confidential information which is proprietary to and a
; trade secret of Visual Numerics Inc.  Use, duplication or disclosure is
; subject to the terms of an appropriate license agreement.
;
   
FUNCTION imsl_csinterp, data1, $              ;INPUT 1-D array: floating point 
                data2, $                   ;INPUT 1-D array: floating point 
                iright=iright, $           ;INPUT Scalar LONG
                right=right, $             ;INPUT Scalar floating point
                ileft=ileft, $             ;INPUT Scalar LONG
                left=left, $               ;INPUT Scalar floating point
                periodic=periodic, $       ;INPUT Scalar ON/OFF flag
                double=double              ;INPUT Scalar ON/OFF flag

@imsl_init.pro
   ON_ERROR, on_err_action
   ; The following checks are made here because 0.0 is a valid value
   ; for the keywords RIGHT and LEFT.

   ;; ITT VIS 
   IF ( KEYWORD_SET(iright))then begin
	if(n_elements(right) EQ 0) THEN $
	   message, 'The keywords IRIGHT and RIGHT must be supplied together.'
   endif
				    
   IF (KEYWORD_SET(ileft))then begin 
       if(n_elements(left) eq 0) THEN $
	     message, 'The keywords ILEFT and LEFT must be supplied together.'
   endif
   nargs = n_params()
   IF (nargs NE 2) THEN message, 'Incorrect number of arguments.'
   RETURN, imsl_compute_spline(data1, data2, $
                          fcn_idx = IMSL_1, $
                          fcn_name="CSINTERP", $
                          iright = iright, $
                          right = right, $
                          ileft = ileft, $
                          left = left, $
                          periodic = periodic, $ 
                          double = double)
END

                   
                   
                   

  
      

  
