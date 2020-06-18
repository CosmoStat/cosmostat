
;--------------------------------------------------------------#
;
; AUTOMATICALLY GENERATED DO NOT MODIFY 

;--------------------------------------------------------------#




;	Wrapping PIOGetTAB3DAxisGrp

FUNCTION PIOGETTAB3DAXISGRP,NAXIS1,NAXIS2,MyGroup
     ON_ERROR,1
     PIOLibIDLSO=shared_lib_path('PIOLibIDL.so')
    if (N_ELEMENTS(NAXIS1) EQ 0) THEN     NAXIS1=LONG64(0) ELSE    NAXIS1=LONG64(NAXIS1)
     if (N_ELEMENTS(NAXIS2) EQ 0) THEN     NAXIS2=LONG64(0) ELSE    NAXIS2=LONG64(NAXIS2)
 
 PIOGETTAB3DAXISGRP=call_external(PIOLibIDLSO,'piogettab3daxisgrp_tempoidl', $
        NAXIS1, $
        NAXIS2, $
        MyGroup, $
               /L64_VALUE) 

 RETURN,PIOGETTAB3DAXISGRP
END
