
;--------------------------------------------------------------#
;
; AUTOMATICALLY GENERATED DO NOT MODIFY 

;--------------------------------------------------------------#




;	Wrapping PIOInfoTAB3DGrp

FUNCTION PIOINFOTAB3DGRP,MyGroup,FLGname=FLGname,NbFLG=NbFLG,TAB3Dtype=TAB3Dtype,TAB3Dname=TAB3Dname,NbTAB3D=NbTAB3D
     ON_ERROR,1
     PIOLibIDLSO=shared_lib_path('PIOLibIDL.so')
    FLGname=long64(0)
    if (N_ELEMENTS(NbFLG) EQ 0) THEN     NbFLG=long(0) ELSE   NbFLG=long(NbFLG)
    TAB3Dtype=long64(0)
    TAB3Dname=long64(0)
    if (N_ELEMENTS(NbTAB3D) EQ 0) THEN     NbTAB3D=long(0) ELSE   NbTAB3D=long(NbTAB3D)

 PIOINFOTAB3DGRP=call_external(PIOLibIDLSO,'pioinfotab3dgrp_tempoidl', $
        FLGname, $
        NbFLG, $
        TAB3Dtype, $
        TAB3Dname, $
        NbTAB3D, $
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
  IF (NbTAB3D GT 0) THEN BEGIN
     TP=TAB3Dtype
     TMPNAME=bytarr(128)
     TAB3Dtype=STRARR(NbTAB3D)
        FOR I=0L,NbTAB3D-1L DO BEGIN
            A=call_external(PIOLibIDLSO,'piogetstringidl',TP+I*128,TMPNAME)
            TAB3Dtype(I)=strmid(string(TMPNAME),0,A)
            END
    END
  IF (NbTAB3D GT 0) THEN BEGIN
     TP=TAB3Dname
     TMPNAME=bytarr(128)
     TAB3Dname=STRARR(NbTAB3D)
        FOR I=0L,NbTAB3D-1L DO BEGIN
            A=call_external(PIOLibIDLSO,'piogetstringidl',TP+I*128,TMPNAME)
            TAB3Dname(I)=strmid(string(TMPNAME),0,A)
            END
    END
    MESSAGE,'*',/CONTINUE
    MESSAGE,'NbTAB3D  :   '+STRING(NbTAB3D),/CONTINUE
    FOR I=0L,NbTAB3D-1L DO MESSAGE,'    '+'   '+STRING(TAB3DTYPE(I))+'   '+STRING(TAB3DNAME(I)) ,/CONTINUE
    MESSAGE,'*',/CONTINUE

 RETURN,PIOINFOTAB3DGRP
END
