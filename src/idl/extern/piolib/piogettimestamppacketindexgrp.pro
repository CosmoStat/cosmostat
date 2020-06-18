
;--------------------------------------------------------------#
;
; AUTOMATICALLY GENERATED DO NOT MODIFY 

;--------------------------------------------------------------#




;	Wrapping PIOGetTIMESTAMPPacketIndexGrp

FUNCTION PIOGETTIMESTAMPPACKETINDEXGRP,Index,MyGroup
     ON_ERROR,1
     PIOLibIDLSO=shared_lib_path('PIOLibIDL.so')
    if (N_ELEMENTS(Index) EQ 0) THEN     Index=LONG64(0) ELSE    Index=LONG64(Index)
 
 PIOGETTIMESTAMPPACKETINDEXGRP=call_external(PIOLibIDLSO,'piogettimestamppacketindexgrp_tempoidl', $
        Index, $
        MyGroup, $
               /L64_VALUE) 

 RETURN,PIOGETTIMESTAMPPACKETINDEXGRP
END
