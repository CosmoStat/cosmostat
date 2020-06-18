; $Id: //depot/idl/IDL_70/idl_src/libs/imsl/imsl_6.0/lib/imsl_cvt_arr.pro#1 $
;
; Copyright (c) 1970-2006, VISUAL NUMERICS Inc. All Rights Reserved.
; This software is confidential information which is proprietary to and a
; trade secret of Visual Numerics Inc.  Use, duplication or disclosure is
; subject to the terms of an appropriate license agreement.
;

FUNCTION imsl_Cvt_arr, Arr, type_desired, back = back
   ; This function assumes that the input array, Arr, is a
   ; valid WAVE 1-D or 2-D array of type float, double,complex, or dcomplex.
   ; In addition to converting the array to the desired data type,
   ; this function also realligns the data in memory for so that
   ; it will be consistent with what C/Math requires.
      
@imsl_init.pro
   ON_ERROR, on_err_action

   ; If BACK is set, then the following cases are handled:
   ;  1. ARR is a complex array:  Just transpose it, and return.
   ;  2. ARR is a double complex array:  Just transpose it, and return.
   ;  3. ARR is a float array: In this case, the data is really complex,
   ;     and we must create a complex array with the correct interleaving
   ;     of the data.
   ;  4. ARR is a double array: In this case, the data is really double 
   ;     complex, and we must create a double complex array with the 
   ;     correct interleaving of the data.
   IF (KEYWORD_SET(back)) THEN BEGIN
      sz_Arr = IMSL_SIZE(Arr)
      dim0 = sz_arr(1)/2
      type_input = IMSL_SIZE(arr, /Type)
      IF (sz_arr(0) EQ 2) THEN BEGIN ; It is a 2-D array...
         dim1 = sz_arr(2)
         CASE type_input OF
            TYP_FLOAT: BEGIN
               ; Otherwise create a single precision complex result
               a_cvt = complexarr(dim1, dim0)
               FOR i = IMSL_0, dim1-1 DO $
                 a_cvt(i, *) = complex(arr(2*lindgen(dim0), i), $
                                       arr(2*lindgen(dim0)+1, i))
               END
            TYP_DOUBLE: BEGIN
               ; If the incoming array is double create a double complex result
               a_cvt = dcomplexarr(dim1, dim0)
               FOR i = IMSL_0, dim1-1 DO $
                 a_cvt(i, *) = dcomplex(arr(2*lindgen(dim0), i), $
                                        arr(2*lindgen(dim0)+1, i))
               END
            TYP_COMPLEX:a_cvt = transpose(Arr)
            TYP_DCMPLX:a_cvt = transpose(Arr)
         END ; Case statement
      END ELSE BEGIN ; It is a 1-D array...
         dim1 = 1
         CASE type_input OF
            TYP_FLOAT: BEGIN
               a_cvt = complexarr(dim0)
               FOR i = IMSL_0, dim1-1 DO $
                 a_cvt(*) = complex(arr(2*lindgen(dim0)), arr(2*lindgen(dim0)+1))
               END
            TYP_DOUBLE: BEGIN
               ; If the incoming array is double create a double complex result
               a_cvt = dcomplexarr(dim0)
               FOR i = IMSL_0, dim1-1 DO $
                 a_cvt(*) = dcomplex(arr(2*lindgen(dim0)), arr(2*lindgen(dim0)+1))
               END
            TYP_COMPLEX:a_cvt = transpose(Arr)
            TYP_DCMPLX:a_cvt = transpose(Arr)
         END ; Case statement
      END
   END ELSE BEGIN
      dim1 = IMSL_1
      sz_Arr = IMSL_SIZE(Arr)
      dim0 = sz_Arr(1)
      IF (sz_Arr(0) EQ 2) THEN dim1 = sz_Arr(2)
      type_Arr = sz_Arr(N_ELEMENTS(sz_Arr)-2)
      CASE type_desired OF
         TYP_FLOAT: a_cvt = float(transpose(Arr))
         TYP_DOUBLE: a_cvt = double(transpose(Arr))
         ;
         ; We have to be careful here since there is the possibility of
         ; trying to convert a dcomplex array to a complex array.
         TYP_COMPLEX:BEGIN
            IF (type_Arr EQ TYP_COMPLEX) THEN a_cvt = transpose(Arr) $
            ELSE BEGIN
               ; Converting a double complex to a complex array.
               IF (type_Arr EQ TYP_DCMPLX) THEN BEGIN
                  a_cvt = fltarr(2*dim1, dim0)
                  FOR i = IMSL_0, dim0-1 DO a_cvt(2*lindgen(dim1), i) = float(Arr(i, *))
                  FOR i = IMSL_0, dim0-1 DO a_cvt(2*lindgen(dim1)+1, i) = imaginary(Arr(i, *))
               END ELSE a_cvt = complex(transpose(Arr)) ; real/double to complex
            END
         END
         TYP_DCMPLX:BEGIN
            IF (type_Arr EQ TYP_DCMPLX) THEN a_cvt = transpose(Arr) $
            ELSE BEGIN 
               ; Converting a  complex to a double complex array.
               IF (type_Arr EQ TYP_COMPLEX) THEN BEGIN
                  a_cvt = dblarr(2*dim1, dim0)
                  FOR i = IMSL_0, dim0-1 DO a_cvt(2*lindgen(dim1), i) = double(Arr(i, *))
                  FOR i = IMSL_0, dim0-1 DO a_cvt(2*lindgen(dim1)+1, i) = imaginary(Arr(i, *))
               END ELSE a_cvt = dcomplex(transpose(Arr))  ; real/double to dcomplex
            END
         END
      END
   END
   RETURN, a_cvt
END

