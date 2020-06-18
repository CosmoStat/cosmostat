
;--------------------------------------------------------------#
;
; AUTOMATICALLY GENERATED DO NOT MODIFY 

;--------------------------------------------------------------#




;	Wrapping PIOInfoTOIGrp

FUNCTION PIOINFOTOIGRP,MyGroup,FLGname=FLGname,NbFLG=NbFLG,TOItype=TOItype,TOIname=TOIname,BeginIndx=BeginIndx,EndIndx=EndIndx,NbTOI=NbTOI,ROIGroup=ROIGroup
     ON_ERROR,1
     PIOLibIDLSO=shared_lib_path('PIOLibIDL.so')
    FLGname=long64(0)
    if (N_ELEMENTS(NbFLG) EQ 0) THEN     NbFLG=long(0) ELSE   NbFLG=long(NbFLG)
    TOItype=long64(0)
    TOIname=long64(0)
    BeginIndx=long64(0)
    EndIndx=long64(0)
    if (N_ELEMENTS(NbTOI) EQ 0) THEN     NbTOI=long(0) ELSE   NbTOI=long(NbTOI)
    ROIGroup_TMP=BYTARR(128)
    ROIGroup_TMP(*)=0
    if (N_ELEMENTS(ROIGroup) GT 0) THEN if (STRLEN(ROIGroup) GT 0) THEN ROIGroup_TMP(0:STRLEN(ROIGroup)-1)=BYTE(ROIGroup)

 PIOINFOTOIGRP=call_external(PIOLibIDLSO,'pioinfotoigrp_tempoidl', $
        FLGname, $
        NbFLG, $
        TOItype, $
        TOIname, $
        BeginIndx, $
        EndIndx, $
        NbTOI, $
        ROIGroup_TMP, $
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
  IF (NbTOI GT 0) THEN BEGIN
     TP=TOItype
     TMPNAME=bytarr(128)
     TOItype=STRARR(NbTOI)
        FOR I=0L,NbTOI-1L DO BEGIN
            A=call_external(PIOLibIDLSO,'piogetstringidl',TP+I*128,TMPNAME)
            TOItype(I)=strmid(string(TMPNAME),0,A)
            END
    END
  IF (NbTOI GT 0) THEN BEGIN
     TP=TOIname
     TMPNAME=bytarr(128)
     TOIname=STRARR(NbTOI)
        FOR I=0L,NbTOI-1L DO BEGIN
            A=call_external(PIOLibIDLSO,'piogetstringidl',TP+I*128,TMPNAME)
            TOIname(I)=strmid(string(TMPNAME),0,A)
            END
    END
  IF (NbTOI GT 0) THEN BEGIN
     TP=BeginIndx
     BeginIndx=LON64ARR(NbTOI)
     A=call_external(PIOLibIDLSO,'piomemcpyidl',TP,BeginIndx, LONG64(NbTOI*8))
    END
  IF (NbTOI GT 0) THEN BEGIN
     TP=EndIndx
     EndIndx=LON64ARR(NbTOI)
     A=call_external(PIOLibIDLSO,'piomemcpyidl',TP,EndIndx, LONG64(NbTOI*8))
    END
    MESSAGE,'*',/CONTINUE
    MESSAGE,'NbTOI  :   '+STRING(NbTOI),/CONTINUE
    FOR I=0L,NbTOI-1L DO MESSAGE,'    '+'   '+STRING(TOITYPE(I))+'   '+STRING(TOINAME(I))+'   '+STRING(BEGININDX(I))+'   '+STRING(ENDINDX(I)) ,/CONTINUE
    MESSAGE,'*',/CONTINUE
  ROIGroup=STRING(ROIGroup_TMP)
    MESSAGE,'*',/CONTINUE
    MESSAGE,' ROIGroup  :'+STRING(ROIGroup) ,/CONTINUE
     MESSAGE,'*',/CONTINUE

 RETURN,PIOINFOTOIGRP
END
