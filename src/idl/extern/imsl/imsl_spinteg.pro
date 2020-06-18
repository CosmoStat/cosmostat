; $Id: //depot/idl/IDL_70/idl_src/libs/imsl/imsl_6.0/lib/imsl_spinteg.pro#1 $   
;
; Copyright (c) 1970-2006, VISUAL NUMERICS Inc. All Rights Reserved.
; This software is confidential information which is proprietary to and a
; trade secret of Visual Numerics Inc.  Use, duplication or disclosure is
; subject to the terms of an appropriate license agreement.
;
FUNCTION imsl_spinteg, data1, $         ;INPUT Scalar floating point 
                    data2, $         ;INPUT Scalar floating point 
                    data3, $         ;INPUT Spline structure, or 
                                     ;      Scalar floating point 
                    data4, $         ;INPUT Scalar floating point 
                    data5            ;INPUT Spline structure,

@imsl_init.pro
   ON_ERROR, on_err_action
   ;
   ; Error checking:
   ; Error checking.
   ;  (argc EQ 3)
   ;    - We are integrating a 1-D spline.
   ;    - DATA1 must be a scalar
   ;    - DATA2 must be a scalar
   ;    - DATA3 must be a spline stucture.
   ;  (argc EQ 5)
   ;    - We are integrating a 2-D spline.
   ;    - DATA1 must be a scalar
   ;    - DATA2 must be a scalar
   ;    - DATA3 must be a scalar
   ;    - DATA4 must be a scalar
   ;    - DATA5 must be a spline stucture.
   ;

   nargs = n_params()
   IF ((nargs NE 3) AND (nargs NE 5)) THEN message, 'Incorrect number of arguments.'
   
   IF (nargs EQ 3) THEN BEGIN
      data3_cvt = imsl_chk_spline_struct(data3, type, struct_is_ppoly)
      domain_dim = IMSL_1
   END ELSE BEGIN
      ; (nargs eq 5)
      data5_cvt = imsl_chk_spline_struct(data5, type, struct_is_ppoly)
      domain_dim = IMSL_2
   END

   ; Error checking complete.
   ;
   ; Decide on what precision to use.
   ; The precision is dependent upon the precision of the data in the
   ; spline structure.  The procedure chk_spline_struct() returns the precision.
   ;
   ; Setup the parameters for the call to the system function.
   ;
   ; Input LONG keyword(s)
   IF (KEYWORD_SET(xderiv)) THEN xderiv_cvt = IMSL_LONG(xderiv(0)) ELSE xderiv_cvt = IMSL_0
   IF (KEYWORD_SET(yderiv)) THEN yderiv_cvt = IMSL_LONG(yderiv(0)) ELSE yderiv_cvt = IMSL_0
   ;
   ; Output LONG keyword(s)
   ;
   ; Floating point arguments and keywords
   IF (type EQ TYP_DOUBLE) THEN BEGIN
      IF (nargs EQ 3) THEN BEGIN
         data1_cvt = double(data1)
         data2_cvt = double(data2)
         result = double(0.)
      END ELSE BEGIN
         data1_cvt = double(data1)
         data2_cvt = double(data2)
         data3_cvt = double(data3)
         data4_cvt = double(data4)
         result =  double(0.)
      END
   END ELSE BEGIN
      IF (nargs EQ 3) THEN BEGIN
         data1_cvt = float(data1)
         data2_cvt = float(data2)
         result = float(0.)
      END ELSE BEGIN 
         data1_cvt = float(data1)
         data2_cvt = float(data2)
         data3_cvt = float(data3)
         data4_cvt = float(data4)
         result =  float(0.)
      END
   END
   ;
   ; Call the system function.
   ;
   err_status = 0L
   MATHSTAT_189, type, err_status, data1_cvt, data2_cvt, data3_cvt, $
                              data4_cvt, data5_cvt, $
                              domain_dim, struct_is_ppoly, $
                              result
                           
   ; return
   RETURN, result
END

                   
                   
                   

  
      

  
