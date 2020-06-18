; $Id: //depot/idl/IDL_70/idl_src/libs/imsl/imsl_6.0/lib/imsl_n_elements.pro#1 $
;
; Copyright (c) 2006-2006, ITT Visual Information Solutions. All
;       rights reserved. Unauthorized reproduction is prohibited.
;
; Helper routine to return N_ELEMENTS as either a
; 32-bit long or 64-bit long, depending on the architecture.
; Required because a PV-Wave type LONG can be either
; 32 or 64 bit.
;
function imsl_n_elements, x

@imsl_init.pro
  ON_ERROR, on_err_action
  return, IMSL64BIT ? LONG64(N_ELEMENTS(x)) : LONG(N_ELEMENTS(x))

end
