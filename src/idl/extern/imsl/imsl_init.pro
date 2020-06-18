; $Id: //depot/idl/IDL_70/idl_src/libs/imsl/imsl_6.0/lib/imsl_init.pro#1 $
;
; Copyright (c) 2006-2006, ITT Visual Information Solutions. All
;       rights reserved. Unauthorized reproduction is prohibited.
;
; IMSL_Init
;
; Purpose:
;    This routine takes the place of the math_init & stat_init routines
;    that ship with wave. This is used to provide an automatic method
;    to init the .pro layer of IMSL. Wave normally requires the user
;    to call math_init and stat_init at the start of a session. This routine
;    unifies these routines and allows for a single entry to be placed at the
;    top of each IMSL routine.
;
; This file should be included with "@imsl_init.pro"
;
    COMMON CMAST$gen_com, on_err_action, FALSE, TRUE, $
        TYP_FLOAT, TYP_DOUBLE, TYP_COMPLEX, TYP_DCMPLX, $
        TYP_LONG, TYP_MEMINT, IMSL64BIT, $
        IMSL_0, IMSL_1, IMSL_2, IMSL_3, IMSL_4

    ; If ON_ERR_ACTION is already defined, just return.
    if (n_elements(ON_ERR_ACTION) eq 0) then begin
        ;Set defaults used in many codes
        ON_ERR_ACTION = 2
        IMSL64BIT = !version.memory_bits eq 64
        ; A PV-Wave type LONG can be either 32 or 64 bit, depending upon
        ; the architecture. Because of this, we need to use the
        ; appropriate constants within the wrappers.
        IMSL_0 = IMSL64BIT ? 0LL : 0L
        IMSL_1 = IMSL64BIT ? 1LL : 1L
        IMSL_2 = IMSL64BIT ? 2LL : 2L
        IMSL_3 = IMSL64BIT ? 3LL : 3L
        IMSL_4 = IMSL64BIT ? 4LL : 4L
        FALSE         = IMSL_0
        TRUE          = IMSL_1
        TYP_LONG      = 3L
        TYP_MEMINT    = IMSL64BIT ? 14L : 3L
        TYP_FLOAT     = 4L
        TYP_DOUBLE    = 5L
        TYP_COMPLEX   = 6L
        TYP_DCMPLX    = 9L
        libPath = FILEPATH('', SUBDIR=['lib','imsl'])
        SETENV, 'IMSLERRPATH=' + libPath
        SETENV, 'IMSLSERRPATH=' + libPath
    endif

; end
