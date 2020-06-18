; $Id: //depot/idl/IDL_70/idl_src/libs/imsl/imsl_6.0/lib/imsl_lonarr.pro#1 $
;
; Copyright (c) 2006-2006, ITT Visual Information Solutions. All
;       rights reserved. Unauthorized reproduction is prohibited.
;
; Helper routine to return a 32-bit or a 64-bit
; long array, depending on the architecture.
; Required because a PV-Wave type LONG can be either
; 32 or 64 bit.
;
function imsl_lonarr, n1, n2, n3

@imsl_init.pro
  ON_ERROR, on_err_action
  case N_Params() of
  0: MESSAGE, 'Incorrect number of arguments.'
  1: return, IMSL64BIT ? LON64ARR(n1) : LONARR(n1)
  2: return, IMSL64BIT ? LON64ARR(n1,n2) : LONARR(n1,n2)
  3: return, IMSL64BIT ? LON64ARR(n1,n2,n3) : LONARR(n1,n2,n3)
  endcase

end
