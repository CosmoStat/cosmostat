
;--------------------------------------------------------------#
;
; AUTOMATICALLY GENERATED DO NOT MODIFY 

;--------------------------------------------------------------#




;	Wrapping PIOWriteRingIndexGrp

FUNCTION PIOWRITERINGINDEXGRP,BeginIndx,EndIndx,RingNum,NbNum,MyGroup
     ON_ERROR,1
     PIOLibIDLSO=shared_lib_path('PIOLibIDL.so')
    if (N_ELEMENTS(BeginIndx) EQ 0) THEN     BeginIndx=LONG64(0) ELSE    BeginIndx=LONG64(BeginIndx)
     if (N_ELEMENTS(EndIndx) EQ 0) THEN     EndIndx=LONG64(0) ELSE    EndIndx=LONG64(EndIndx)
     if (N_ELEMENTS(RingNum) EQ 0) THEN     RingNum=LONG64(0) ELSE    RingNum=LONG64(RingNum)
     if (N_ELEMENTS(NbNum) EQ 0) THEN     NbNum=LONG64(0) ELSE    NbNum=LONG64(NbNum)
 
MEM_TYP=[-1,1,2,4,4,8,8,-1,-1,16,-1,-1,-1,-1,8,-1]
 PIOWRITERINGINDEXGRP=call_external(PIOLibIDLSO,'piowriteringindexgrp_tempoidl', $
        BeginIndx, $
        LONG64(N_ELEMENTS(BeginIndx)*MEM_TYP(size(BeginIndx,/type))), $
        EndIndx, $
        LONG64(N_ELEMENTS(EndIndx)*MEM_TYP(size(EndIndx,/type))), $
        RingNum, $
        LONG64(N_ELEMENTS(RingNum)*MEM_TYP(size(RingNum,/type))), $
        NbNum, $
        MyGroup, $
               /L64_VALUE) 

 RETURN,PIOWRITERINGINDEXGRP
END
