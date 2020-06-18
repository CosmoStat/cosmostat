
;--------------------------------------------------------------#
;
; AUTOMATICALLY GENERATED DO NOT MODIFY 

;--------------------------------------------------------------#




;	Wrapping PIOInfoTIMESTAMPGrp

FUNCTION PIOINFOTIMESTAMPGRP,MyGroup,TIMESTAMPtype=TIMESTAMPtype,TIMESTAMPname=TIMESTAMPname,BeginIndx=BeginIndx,EndIndx=EndIndx,NbTIMESTAMP=NbTIMESTAMP,PacketGroup=PacketGroup
     ON_ERROR,1
     PIOLibIDLSO=shared_lib_path('PIOLibIDL.so')
    TIMESTAMPtype=long64(0)
    TIMESTAMPname=long64(0)
    BeginIndx=long64(0)
    EndIndx=long64(0)
    if (N_ELEMENTS(NbTIMESTAMP) EQ 0) THEN     NbTIMESTAMP=long(0) ELSE   NbTIMESTAMP=long(NbTIMESTAMP)
    PacketGroup_TMP=BYTARR(128)
    PacketGroup_TMP(*)=0
    if (N_ELEMENTS(PacketGroup) GT 0) THEN if (STRLEN(PacketGroup) GT 0) THEN PacketGroup_TMP(0:STRLEN(PacketGroup)-1)=BYTE(PacketGroup)

 PIOINFOTIMESTAMPGRP=call_external(PIOLibIDLSO,'pioinfotimestampgrp_tempoidl', $
        TIMESTAMPtype, $
        TIMESTAMPname, $
        BeginIndx, $
        EndIndx, $
        NbTIMESTAMP, $
        PacketGroup_TMP, $
        MyGroup, $
               /L64_VALUE) 
  IF (NbTIMESTAMP GT 0) THEN BEGIN
     TP=TIMESTAMPtype
     TMPNAME=bytarr(128)
     TIMESTAMPtype=STRARR(NbTIMESTAMP)
        FOR I=0L,NbTIMESTAMP-1L DO BEGIN
            A=call_external(PIOLibIDLSO,'piogetstringidl',TP+I*128,TMPNAME)
            TIMESTAMPtype(I)=strmid(string(TMPNAME),0,A)
            END
    END
  IF (NbTIMESTAMP GT 0) THEN BEGIN
     TP=TIMESTAMPname
     TMPNAME=bytarr(128)
     TIMESTAMPname=STRARR(NbTIMESTAMP)
        FOR I=0L,NbTIMESTAMP-1L DO BEGIN
            A=call_external(PIOLibIDLSO,'piogetstringidl',TP+I*128,TMPNAME)
            TIMESTAMPname(I)=strmid(string(TMPNAME),0,A)
            END
    END
  IF (NbTIMESTAMP GT 0) THEN BEGIN
     TP=BeginIndx
     BeginIndx=LON64ARR(NbTIMESTAMP)
     A=call_external(PIOLibIDLSO,'piomemcpyidl',TP,BeginIndx, LONG64(NbTIMESTAMP*8))
    END
  IF (NbTIMESTAMP GT 0) THEN BEGIN
     TP=EndIndx
     EndIndx=LON64ARR(NbTIMESTAMP)
     A=call_external(PIOLibIDLSO,'piomemcpyidl',TP,EndIndx, LONG64(NbTIMESTAMP*8))
    END
    MESSAGE,'*',/CONTINUE
    MESSAGE,'NbTIMESTAMP  :   '+STRING(NbTIMESTAMP),/CONTINUE
    FOR I=0L,NbTIMESTAMP-1L DO MESSAGE,'    '+'   '+STRING(TIMESTAMPTYPE(I))+'   '+STRING(TIMESTAMPNAME(I))+'   '+STRING(BEGININDX(I))+'   '+STRING(ENDINDX(I)) ,/CONTINUE
    MESSAGE,'*',/CONTINUE
  PacketGroup=STRING(PacketGroup_TMP)
    MESSAGE,'*',/CONTINUE
    MESSAGE,' PacketGroup  :'+STRING(PacketGroup) ,/CONTINUE
     MESSAGE,'*',/CONTINUE

 RETURN,PIOINFOTIMESTAMPGRP
END
