
;--------------------------------------------------------------#
;
; AUTOMATICALLY GENERATED DO NOT MODIFY 

;--------------------------------------------------------------#




;	Wrapping PIOFreeInfoTOIGrp

FUNCTION PIOFREEINFOTOIGRP,FLGname,TOItype,TOIname,BeginIndx,EndIndx
     ON_ERROR,1
     PIOLibIDLSO=shared_lib_path('PIOLibIDL.so')
    if (N_ELEMENTS(BeginIndx) EQ 0) THEN     BeginIndx=LONG64(0) ELSE    BeginIndx=LONG64(BeginIndx)
     if (N_ELEMENTS(EndIndx) EQ 0) THEN     EndIndx=LONG64(0) ELSE    EndIndx=LONG64(EndIndx)
 
 PIOFREEINFOTOIGRP=call_external(PIOLibIDLSO,'piofreeinfotoigrp_tempoidl', $
        FLGname, $
        TOItype, $
        TOIname, $
        BeginIndx, $
        EndIndx, $
               /L64_VALUE) 

 RETURN,PIOFREEINFOTOIGRP
END
