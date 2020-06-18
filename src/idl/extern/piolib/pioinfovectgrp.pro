
;--------------------------------------------------------------#
;
; AUTOMATICALLY GENERATED DO NOT MODIFY 

;--------------------------------------------------------------#




;	Wrapping PIOInfoVECTGrp

FUNCTION PIOINFOVECTGRP,MyGroup,FLGname=FLGname,NbFLG=NbFLG,VECTtype=VECTtype,VECTname=VECTname,NbVECT=NbVECT
     ON_ERROR,1
     PIOLibIDLSO=shared_lib_path('PIOLibIDL.so')
    FLGname=long64(0)
    if (N_ELEMENTS(NbFLG) EQ 0) THEN     NbFLG=long(0) ELSE   NbFLG=long(NbFLG)
    VECTtype=long64(0)
    VECTname=long64(0)
    if (N_ELEMENTS(NbVECT) EQ 0) THEN     NbVECT=long(0) ELSE   NbVECT=long(NbVECT)

 PIOINFOVECTGRP=call_external(PIOLibIDLSO,'pioinfovectgrp_tempoidl', $
        FLGname, $
        NbFLG, $
        VECTtype, $
        VECTname, $
        NbVECT, $
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
  IF (NbVECT GT 0) THEN BEGIN
     TP=VECTtype
     TMPNAME=bytarr(128)
     VECTtype=STRARR(NbVECT)
        FOR I=0L,NbVECT-1L DO BEGIN
            A=call_external(PIOLibIDLSO,'piogetstringidl',TP+I*128,TMPNAME)
            VECTtype(I)=strmid(string(TMPNAME),0,A)
            END
    END
  IF (NbVECT GT 0) THEN BEGIN
     TP=VECTname
     TMPNAME=bytarr(128)
     VECTname=STRARR(NbVECT)
        FOR I=0L,NbVECT-1L DO BEGIN
            A=call_external(PIOLibIDLSO,'piogetstringidl',TP+I*128,TMPNAME)
            VECTname(I)=strmid(string(TMPNAME),0,A)
            END
    END
    MESSAGE,'*',/CONTINUE
    MESSAGE,'NbVECT  :   '+STRING(NbVECT),/CONTINUE
    FOR I=0L,NbVECT-1L DO MESSAGE,'    '+'   '+STRING(VECTTYPE(I))+'   '+STRING(VECTNAME(I)) ,/CONTINUE
    MESSAGE,'*',/CONTINUE

 RETURN,PIOINFOVECTGRP
END
