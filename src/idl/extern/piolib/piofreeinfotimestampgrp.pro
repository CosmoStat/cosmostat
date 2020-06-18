
;--------------------------------------------------------------#
;
; AUTOMATICALLY GENERATED DO NOT MODIFY 

;--------------------------------------------------------------#




;	Wrapping PIOFreeInfoTIMESTAMPGrp

FUNCTION PIOFREEINFOTIMESTAMPGRP,TIMESTAMPtype,TIMESTAMPname,BeginIndx,EndIndx
     ON_ERROR,1
     PIOLibIDLSO=shared_lib_path('PIOLibIDL.so')
    if (N_ELEMENTS(BeginIndx) EQ 0) THEN     BeginIndx=LONG64(0) ELSE    BeginIndx=LONG64(BeginIndx)
     if (N_ELEMENTS(EndIndx) EQ 0) THEN     EndIndx=LONG64(0) ELSE    EndIndx=LONG64(EndIndx)
 
 PIOFREEINFOTIMESTAMPGRP=call_external(PIOLibIDLSO,'piofreeinfotimestampgrp_tempoidl', $
        TIMESTAMPtype, $
        TIMESTAMPname, $
        BeginIndx, $
        EndIndx, $
               /L64_VALUE) 

 RETURN,PIOFREEINFOTIMESTAMPGRP
END
