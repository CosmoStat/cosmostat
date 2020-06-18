
;--------------------------------------------------------------#
;
; AUTOMATICALLY GENERATED DO NOT MODIFY 

;--------------------------------------------------------------#




;	Wrapping PIOInfoPBRGrp

FUNCTION PIOINFOPBRGRP,MyGroup,FLGname=FLGname,NbFLG=NbFLG,POItype=POItype,POIname=POIname,BeginIndx=BeginIndx,EndIndx=EndIndx,NbPOI=NbPOI,PHItype=PHItype,PHIname=PHIname,NbPHI=NbPHI,ROIGroup=ROIGroup,RingSize=RingSize
     ON_ERROR,1
     PIOLibIDLSO=shared_lib_path('PIOLibIDL.so')
    FLGname=long64(0)
    if (N_ELEMENTS(NbFLG) EQ 0) THEN     NbFLG=long(0) ELSE   NbFLG=long(NbFLG)
    POItype=long64(0)
    POIname=long64(0)
    BeginIndx=long64(0)
    EndIndx=long64(0)
    if (N_ELEMENTS(NbPOI) EQ 0) THEN     NbPOI=long(0) ELSE   NbPOI=long(NbPOI)
    PHItype=long64(0)
    PHIname=long64(0)
    if (N_ELEMENTS(NbPHI) EQ 0) THEN     NbPHI=long(0) ELSE   NbPHI=long(NbPHI)
    ROIGroup_TMP=BYTARR(128)
    ROIGroup_TMP(*)=0
    if (N_ELEMENTS(ROIGroup) GT 0) THEN if (STRLEN(ROIGroup) GT 0) THEN ROIGroup_TMP(0:STRLEN(ROIGroup)-1)=BYTE(ROIGroup)
    if (N_ELEMENTS(RingSize) EQ 0) THEN     RingSize=LONG64(0) ELSE    RingSize=LONG64(RingSize)
 
 PIOINFOPBRGRP=call_external(PIOLibIDLSO,'pioinfopbrgrp_tempoidl', $
        FLGname, $
        NbFLG, $
        POItype, $
        POIname, $
        BeginIndx, $
        EndIndx, $
        NbPOI, $
        PHItype, $
        PHIname, $
        NbPHI, $
        ROIGroup_TMP, $
        RingSize, $
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
  IF (NbPOI GT 0) THEN BEGIN
     TP=POItype
     TMPNAME=bytarr(128)
     POItype=STRARR(NbPOI)
        FOR I=0L,NbPOI-1L DO BEGIN
            A=call_external(PIOLibIDLSO,'piogetstringidl',TP+I*128,TMPNAME)
            POItype(I)=strmid(string(TMPNAME),0,A)
            END
    END
  IF (NbPOI GT 0) THEN BEGIN
     TP=POIname
     TMPNAME=bytarr(128)
     POIname=STRARR(NbPOI)
        FOR I=0L,NbPOI-1L DO BEGIN
            A=call_external(PIOLibIDLSO,'piogetstringidl',TP+I*128,TMPNAME)
            POIname(I)=strmid(string(TMPNAME),0,A)
            END
    END
  IF (NbPOI GT 0) THEN BEGIN
     TP=BeginIndx
     BeginIndx=LON64ARR(NbPOI)
     A=call_external(PIOLibIDLSO,'piomemcpyidl',TP,BeginIndx, LONG64(NbPOI*8))
    END
  IF (NbPOI GT 0) THEN BEGIN
     TP=EndIndx
     EndIndx=LON64ARR(NbPOI)
     A=call_external(PIOLibIDLSO,'piomemcpyidl',TP,EndIndx, LONG64(NbPOI*8))
    END
    MESSAGE,'*',/CONTINUE
    MESSAGE,'NbPOI  :   '+STRING(NbPOI),/CONTINUE
    FOR I=0L,NbPOI-1L DO MESSAGE,'    '+'   '+STRING(POITYPE(I))+'   '+STRING(POINAME(I))+'   '+STRING(BEGININDX(I))+'   '+STRING(ENDINDX(I)) ,/CONTINUE
    MESSAGE,'*',/CONTINUE
  IF (NbPHI GT 0) THEN BEGIN
     TP=PHItype
     TMPNAME=bytarr(128)
     PHItype=STRARR(NbPHI)
        FOR I=0L,NbPHI-1L DO BEGIN
            A=call_external(PIOLibIDLSO,'piogetstringidl',TP+I*128,TMPNAME)
            PHItype(I)=strmid(string(TMPNAME),0,A)
            END
    END
  IF (NbPHI GT 0) THEN BEGIN
     TP=PHIname
     TMPNAME=bytarr(128)
     PHIname=STRARR(NbPHI)
        FOR I=0L,NbPHI-1L DO BEGIN
            A=call_external(PIOLibIDLSO,'piogetstringidl',TP+I*128,TMPNAME)
            PHIname(I)=strmid(string(TMPNAME),0,A)
            END
    END
    MESSAGE,'*',/CONTINUE
    MESSAGE,'NbPHI  :   '+STRING(NbPHI),/CONTINUE
    FOR I=0L,NbPHI-1L DO MESSAGE,'    '+'   '+STRING(PHITYPE(I))+'   '+STRING(PHINAME(I)) ,/CONTINUE
    MESSAGE,'*',/CONTINUE
  ROIGroup=STRING(ROIGroup_TMP)
    MESSAGE,'*',/CONTINUE
    MESSAGE,' ROIGroup  :'+STRING(ROIGroup) ,/CONTINUE
     MESSAGE,'*',/CONTINUE

 RETURN,PIOINFOPBRGRP
END
