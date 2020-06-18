; Helper routine used by CONVOL1D and CORR1D.
;
; Copyright (c) 1970-2006, VISUAL NUMERICS Inc. All Rights Reserved.
; This software is confidential information which is proprietary to and a
; trade secret of Visual Numerics Inc.  Use, duplication or disclosure is
; subject to the terms of an appropriate license agreement.
;
;
FUNCTION imsl_pow235, nx, ny, ipad
@imsl_init.pro
  ON_ERROR, on_err_action
   ;
   ; $Id: //depot/idl/IDL_70/idl_src/libs/imsl/imsl_6.0/lib/imsl_pow235.pro#1 $
   ;
   ; If ipad is zero, then return MAX(nx, ny), otherwise return the 
   ; smallest long integer of the form n = (2^a)*(3^b)*(5^c) for nonzero 
   ; whole numbers a, b, c such that n GE (nx+ny-1)
   ;
   IF NOT KEYWORD_SET(ny) THEN ny = nx
   IF (ipad EQ 0) THEN n = IMSL_LONG(nx > ny) $
     ELSE BEGIN 
      nz2 = double(nx+ny-1)
      valold = 2.*nz2
      iend =  IMSL_LONG(alog(nz2)/alog(5.0e0) + .5)
      FOR i = 0, iend DO  BEGIN 
         jend =  IMSL_LONG(alog(nz2/(5.^i))/alog(3.) + .5)
         FOR j = 0, jend DO  BEGIN 
            kstart = IMSL_LONG(alog(nz2/((5.^i)*(3.^j)))/alog(2.) + .5)
            kend = IMSL_LONG((alog((float(nz2)/(5.0e0^i))/(3.0e0^j)))/ $
                          (alog(2.0e0)) + .5)
            FOR k = kstart, kend DO BEGIN
               valnew = (2.^k)*(3.^j)*(5.^i)
               IF ((valnew GE nz2) AND (valnew LT valold)) THEN BEGIN 
                 valold = valnew
                 IF (valold EQ  nz2) THEN goto, FOUND
              END ;if statement
           END ; k loop
        END ; j loop
     END ; i loop
  END ; if (ipad ne 0)...
 
  FOUND:   IF (ipad EQ 1) THEN n = valold
  RETURN, IMSL_LONG(n)
END
