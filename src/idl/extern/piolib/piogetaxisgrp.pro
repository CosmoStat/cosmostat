
;--------------------------------------------------------------#
;
; AUTOMATICALLY GENERATED DO NOT MODIFY 

;--------------------------------------------------------------#




;	Wrapping PIOGetAxisGrp

FUNCTION PIOGETAXISGRP,NAXIS1,NAXIS2,NAXIS3,MyGroup
     ON_ERROR,1
     PIOLibIDLSO=shared_lib_path('PIOLibIDL.so')
    if (N_ELEMENTS(NAXIS1) EQ 0) THEN     NAXIS1=LONG64(0) ELSE    NAXIS1=LONG64(NAXIS1)
     if (N_ELEMENTS(NAXIS2) EQ 0) THEN     NAXIS2=LONG64(0) ELSE    NAXIS2=LONG64(NAXIS2)
     if (N_ELEMENTS(NAXIS3) EQ 0) THEN     NAXIS3=LONG64(0) ELSE    NAXIS3=LONG64(NAXIS3)
 
 PIOGETAXISGRP=call_external(PIOLibIDLSO,'piogetaxisgrp_tempoidl', $
        NAXIS1, $
        NAXIS2, $
        NAXIS3, $
        MyGroup, $
               /L64_VALUE) 

 RETURN,PIOGETAXISGRP
END
