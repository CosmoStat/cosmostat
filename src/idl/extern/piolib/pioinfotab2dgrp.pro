
;--------------------------------------------------------------#
;
; AUTOMATICALLY GENERATED DO NOT MODIFY 

;--------------------------------------------------------------#




;	Wrapping PIOInfoTAB2DGrp

FUNCTION PIOINFOTAB2DGRP,MyGroup,FLGname=FLGname,NbFLG=NbFLG,TAB2Dtype=TAB2Dtype,TAB2Dname=TAB2Dname,NbTAB2D=NbTAB2D
     ON_ERROR,1
     PIOLibIDLSO=shared_lib_path('PIOLibIDL.so')
    FLGname=long64(0)
    if (N_ELEMENTS(NbFLG) EQ 0) THEN     NbFLG=long(0) ELSE   NbFLG=long(NbFLG)
    TAB2Dtype=long64(0)
    TAB2Dname=long64(0)
    if (N_ELEMENTS(NbTAB2D) EQ 0) THEN     NbTAB2D=long(0) ELSE   NbTAB2D=long(NbTAB2D)

 PIOINFOTAB2DGRP=call_external(PIOLibIDLSO,'pioinfotab2dgrp_tempoidl', $
        FLGname, $
        NbFLG, $
        TAB2Dtype, $
        TAB2Dname, $
        NbTAB2D, $
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
  IF (NbTAB2D GT 0) THEN BEGIN
     TP=TAB2Dtype
     TMPNAME=bytarr(128)
     TAB2Dtype=STRARR(NbTAB2D)
        FOR I=0L,NbTAB2D-1L DO BEGIN
            A=call_external(PIOLibIDLSO,'piogetstringidl',TP+I*128,TMPNAME)
            TAB2Dtype(I)=strmid(string(TMPNAME),0,A)
            END
    END
  IF (NbTAB2D GT 0) THEN BEGIN
     TP=TAB2Dname
     TMPNAME=bytarr(128)
     TAB2Dname=STRARR(NbTAB2D)
        FOR I=0L,NbTAB2D-1L DO BEGIN
            A=call_external(PIOLibIDLSO,'piogetstringidl',TP+I*128,TMPNAME)
            TAB2Dname(I)=strmid(string(TMPNAME),0,A)
            END
    END
    MESSAGE,'*',/CONTINUE
    MESSAGE,'NbTAB2D  :   '+STRING(NbTAB2D),/CONTINUE
    FOR I=0L,NbTAB2D-1L DO MESSAGE,'    '+'   '+STRING(TAB2DTYPE(I))+'   '+STRING(TAB2DNAME(I)) ,/CONTINUE
    MESSAGE,'*',/CONTINUE

 RETURN,PIOINFOTAB2DGRP
END
