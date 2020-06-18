
;--------------------------------------------------------------#
;
; AUTOMATICALLY GENERATED DO NOT MODIFY 

;--------------------------------------------------------------#




;	Wrapping PIOGetPVGrp

FUNCTION PIOGETPVGRP,PV,MyGroup
     ON_ERROR,1
     PIOLibIDLSO=shared_lib_path('PIOLibIDL.so')
    if (N_ELEMENTS(PV) EQ 0) THEN     if (N_ELEMENTS(PV) EQ 0) THEN     if (N_ELEMENTS(PV) EQ 0) THEN 
 PIOGETPVGRP=call_external(PIOLibIDLSO,'piogetpvgrp_tempoidl', $
        PV, $
        MyGroup, $
               /L64_VALUE) 

 RETURN,PIOGETPVGRP
END
