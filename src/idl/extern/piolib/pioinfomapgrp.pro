
;--------------------------------------------------------------#
;
; AUTOMATICALLY GENERATED DO NOT MODIFY 

;--------------------------------------------------------------#




;	Wrapping PIOInfoMAPGrp

FUNCTION PIOINFOMAPGRP,MyGroup,FLGname=FLGname,NbFLG=NbFLG,MAPtype=MAPtype,MAPname=MAPname,NbMAP=NbMAP,ordering=ordering,coordsys=coordsys,Nside=Nside
     ON_ERROR,1
     PIOLibIDLSO=shared_lib_path('PIOLibIDL.so')
    FLGname=long64(0)
    if (N_ELEMENTS(NbFLG) EQ 0) THEN     NbFLG=long(0) ELSE   NbFLG=long(NbFLG)
    MAPtype=long64(0)
    MAPname=long64(0)
    if (N_ELEMENTS(NbMAP) EQ 0) THEN     NbMAP=long(0) ELSE   NbMAP=long(NbMAP)
    ordering_TMP=BYTARR(128)
    ordering_TMP(*)=0
    if (N_ELEMENTS(ordering) GT 0) THEN if (STRLEN(ordering) GT 0) THEN ordering_TMP(0:STRLEN(ordering)-1)=BYTE(ordering)
    coordsys_TMP=BYTARR(128)
    coordsys_TMP(*)=0
    if (N_ELEMENTS(coordsys) GT 0) THEN if (STRLEN(coordsys) GT 0) THEN coordsys_TMP(0:STRLEN(coordsys)-1)=BYTE(coordsys)
    if (N_ELEMENTS(Nside) EQ 0) THEN     Nside=long(0) ELSE   Nside=long(Nside)

 PIOINFOMAPGRP=call_external(PIOLibIDLSO,'pioinfomapgrp_tempoidl', $
        FLGname, $
        NbFLG, $
        MAPtype, $
        MAPname, $
        NbMAP, $
        ordering_TMP, $
        coordsys_TMP, $
        Nside, $
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
  IF (NbMAP GT 0) THEN BEGIN
     TP=MAPtype
     TMPNAME=bytarr(128)
     MAPtype=STRARR(NbMAP)
        FOR I=0L,NbMAP-1L DO BEGIN
            A=call_external(PIOLibIDLSO,'piogetstringidl',TP+I*128,TMPNAME)
            MAPtype(I)=strmid(string(TMPNAME),0,A)
            END
    END
  IF (NbMAP GT 0) THEN BEGIN
     TP=MAPname
     TMPNAME=bytarr(128)
     MAPname=STRARR(NbMAP)
        FOR I=0L,NbMAP-1L DO BEGIN
            A=call_external(PIOLibIDLSO,'piogetstringidl',TP+I*128,TMPNAME)
            MAPname(I)=strmid(string(TMPNAME),0,A)
            END
    END
    MESSAGE,'*',/CONTINUE
    MESSAGE,'NbMAP  :   '+STRING(NbMAP),/CONTINUE
    FOR I=0L,NbMAP-1L DO MESSAGE,'    '+'   '+STRING(MAPTYPE(I))+'   '+STRING(MAPNAME(I)) ,/CONTINUE
    MESSAGE,'*',/CONTINUE
  ordering=STRING(ordering_TMP)
    MESSAGE,'*',/CONTINUE
    MESSAGE,' ordering  :'+STRING(ordering) ,/CONTINUE
     MESSAGE,'*',/CONTINUE
  coordsys=STRING(coordsys_TMP)
    MESSAGE,'*',/CONTINUE
    MESSAGE,' coordsys  :'+STRING(coordsys) ,/CONTINUE
     MESSAGE,'*',/CONTINUE

 RETURN,PIOINFOMAPGRP
END
