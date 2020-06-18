
;--------------------------------------------------------------#
;
; AUTOMATICALLY GENERATED DO NOT MODIFY 

;--------------------------------------------------------------#




;	Wrapping PIOReadRingIndexGrp

FUNCTION PIOREADRINGINDEXGRP,BeginIndx,EndIndx,RingNum,MyGroup
     ON_ERROR,1
     PIOLibIDLSO=shared_lib_path('PIOLibIDL.so')
    BeginIndx=long64(0)
    EndIndx=long64(0)
    RingNum=long64(0)

 PIOREADRINGINDEXGRP=call_external(PIOLibIDLSO,'pioreadringindexgrp_tempoidl', $
        BeginIndx, $
        EndIndx, $
        RingNum, $
        MyGroup, $
               /L64_VALUE) 
  IF (PIOREADRINGINDEXGRP GT 0) THEN BEGIN
     TP=BeginIndx
     BeginIndx=LON64ARR(PIOREADRINGINDEXGRP)
     A=call_external(PIOLibIDLSO,'piomemcpyidl',TP,BeginIndx, LONG64(PIOREADRINGINDEXGRP*8))
    END
  IF (PIOREADRINGINDEXGRP GT 0) THEN BEGIN
     TP=EndIndx
     EndIndx=LON64ARR(PIOREADRINGINDEXGRP)
     A=call_external(PIOLibIDLSO,'piomemcpyidl',TP,EndIndx, LONG64(PIOREADRINGINDEXGRP*8))
    END
  IF (PIOREADRINGINDEXGRP GT 0) THEN BEGIN
     TP=RingNum
     RingNum=LON64ARR(PIOREADRINGINDEXGRP)
     A=call_external(PIOLibIDLSO,'piomemcpyidl',TP,RingNum, LONG64(PIOREADRINGINDEXGRP*8))
    END

 RETURN,PIOREADRINGINDEXGRP
END
