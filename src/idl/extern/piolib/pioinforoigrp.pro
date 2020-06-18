
;--------------------------------------------------------------#
;
; AUTOMATICALLY GENERATED DO NOT MODIFY 

;--------------------------------------------------------------#




;	Wrapping PIOInfoROIGrp

FUNCTION PIOINFOROIGRP,MyGroup,FLGname=FLGname,NbFLG=NbFLG,ROItype=ROItype,ROIname=ROIname,BeginIndx=BeginIndx,EndIndx=EndIndx,NbROI=NbROI
     ON_ERROR,1
     PIOLibIDLSO=shared_lib_path('PIOLibIDL.so')
    FLGname=long64(0)
    if (N_ELEMENTS(NbFLG) EQ 0) THEN     NbFLG=long(0) ELSE   NbFLG=long(NbFLG)
    ROItype=long64(0)
    ROIname=long64(0)
    BeginIndx=long64(0)
    EndIndx=long64(0)
    if (N_ELEMENTS(NbROI) EQ 0) THEN     NbROI=long(0) ELSE   NbROI=long(NbROI)

 PIOINFOROIGRP=call_external(PIOLibIDLSO,'pioinforoigrp_tempoidl', $
        FLGname, $
        NbFLG, $
        ROItype, $
        ROIname, $
        BeginIndx, $
        EndIndx, $
        NbROI, $
        MyGroup, $
               /L64_VALUE) 
  IF (NbFLG GT 0) THEN BEGIN
     TP=FLGname
     TMPNAME=bytarr(128)
     FLGname=STRARR(NbFLG)
        FOR I=0L,NbFLG-1L DO BEGIN
            A=call_external(PIOLibIDLSO,'piogetstringidl',TP+I*128,TMPNAME)
            FLGname(I)=strmid(string(TMPNAME),0,A)
            END
    END
    MESSAGE,'*',/CONTINUE
    MESSAGE,'NbFLG  :   '+STRING(NbFLG),/CONTINUE
    FOR I=0L,NbFLG-1L DO MESSAGE,'    '+'   '+STRING(FLGNAME(I)) ,/CONTINUE
    MESSAGE,'*',/CONTINUE
  IF (NbROI GT 0) THEN BEGIN
     TP=ROItype
     TMPNAME=bytarr(128)
     ROItype=STRARR(NbROI)
        FOR I=0L,NbROI-1L DO BEGIN
            A=call_external(PIOLibIDLSO,'piogetstringidl',TP+I*128,TMPNAME)
            ROItype(I)=strmid(string(TMPNAME),0,A)
            END
    END
  IF (NbROI GT 0) THEN BEGIN
     TP=ROIname
     TMPNAME=bytarr(128)
     ROIname=STRARR(NbROI)
        FOR I=0L,NbROI-1L DO BEGIN
            A=call_external(PIOLibIDLSO,'piogetstringidl',TP+I*128,TMPNAME)
            ROIname(I)=strmid(string(TMPNAME),0,A)
            END
    END
  IF (NbROI GT 0) THEN BEGIN
     TP=BeginIndx
     BeginIndx=LON64ARR(NbROI)
     A=call_external(PIOLibIDLSO,'piomemcpyidl',TP,BeginIndx, LONG64(NbROI*8))
    END
  IF (NbROI GT 0) THEN BEGIN
     TP=EndIndx
     EndIndx=LON64ARR(NbROI)
     A=call_external(PIOLibIDLSO,'piomemcpyidl',TP,EndIndx, LONG64(NbROI*8))
    END
    MESSAGE,'*',/CONTINUE
    MESSAGE,'NbROI  :   '+STRING(NbROI),/CONTINUE
    FOR I=0L,NbROI-1L DO MESSAGE,'    '+'   '+STRING(ROITYPE(I))+'   '+STRING(ROINAME(I))+'   '+STRING(BEGININDX(I))+'   '+STRING(ENDINDX(I)) ,/CONTINUE
    MESSAGE,'*',/CONTINUE

 RETURN,PIOINFOROIGRP
END
