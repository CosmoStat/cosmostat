
;--------------------------------------------------------------#
;
; AUTOMATICALLY GENERATED DO NOT MODIFY 

;--------------------------------------------------------------#




;	Wrapping PIOGetTOIRingIndexGrp

FUNCTION PIOGETTOIRINGINDEXGRP,Index,MyGroup
     ON_ERROR,1
     PIOLibIDLSO=shared_lib_path('PIOLibIDL.so')
    if (N_ELEMENTS(Index) EQ 0) THEN     Index=LONG64(0) ELSE    Index=LONG64(Index)
 
 PIOGETTOIRINGINDEXGRP=call_external(PIOLibIDLSO,'piogettoiringindexgrp_tempoidl', $
        Index, $
        MyGroup, $
               /L64_VALUE) 

 RETURN,PIOGETTOIRINGINDEXGRP
END
