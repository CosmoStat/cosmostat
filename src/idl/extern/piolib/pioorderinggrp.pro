
;--------------------------------------------------------------#
;
; AUTOMATICALLY GENERATED DO NOT MODIFY 

;--------------------------------------------------------------#




;	Wrapping PIOOrderingGrp

FUNCTION PIOORDERINGGRP,ordering,MyGroup
     ON_ERROR,1
     PIOLibIDLSO=shared_lib_path('PIOLibIDL.so')
    ordering_TMP=BYTARR(128)
    ordering_TMP(*)=0
    if (N_ELEMENTS(ordering) GT 0) THEN if (STRLEN(ordering) GT 0) THEN ordering_TMP(0:STRLEN(ordering)-1)=BYTE(ordering)

 PIOORDERINGGRP=call_external(PIOLibIDLSO,'pioorderinggrp_tempoidl', $
        ordering_TMP, $
        MyGroup, $
               /L64_VALUE) 
  ordering=STRING(ordering_TMP)

 RETURN,PIOORDERINGGRP
END
