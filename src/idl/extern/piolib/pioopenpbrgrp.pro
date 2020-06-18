
;--------------------------------------------------------------#
;
; AUTOMATICALLY GENERATED DO NOT MODIFY 

;--------------------------------------------------------------#




;	Wrapping PIOOpenPBRGrp

FUNCTION PIOOPENPBRGRP,object,mode
     ON_ERROR,1
     PIOLibIDLSO=shared_lib_path('PIOLibIDL.so')
    object_TMP=BYTARR(128)
    object_TMP(*)=0
    if (N_ELEMENTS(object) GT 0) THEN if (STRLEN(object) GT 0) THEN object_TMP(0:STRLEN(object)-1)=BYTE(object)
    mode_TMP=BYTARR(128)
    mode_TMP(*)=0
    if (N_ELEMENTS(mode) GT 0) THEN if (STRLEN(mode) GT 0) THEN mode_TMP(0:STRLEN(mode)-1)=BYTE(mode)

 PIOOPENPBRGRP=call_external(PIOLibIDLSO,'pioopenpbrgrp_tempoidl', $
        object_TMP, $
        mode_TMP, $
               /L64_VALUE) 

 RETURN,PIOOPENPBRGRP
END
