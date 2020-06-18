
;--------------------------------------------------------------#
;
; AUTOMATICALLY GENERATED DO NOT MODIFY 

;--------------------------------------------------------------#




;	Wrapping PIOInfoIMG2DGrp

FUNCTION PIOINFOIMG2DGRP,MyGroup,FLGname=FLGname,NbFLG=NbFLG,IMG2Dtype=IMG2Dtype,IMG2Dname=IMG2Dname,NbIMG2D=NbIMG2D
     ON_ERROR,1
     PIOLibIDLSO=shared_lib_path('PIOLibIDL.so')
    FLGname=long64(0)
    if (N_ELEMENTS(NbFLG) EQ 0) THEN     NbFLG=long(0) ELSE   NbFLG=long(NbFLG)
    IMG2Dtype=long64(0)
    IMG2Dname=long64(0)
    if (N_ELEMENTS(NbIMG2D) EQ 0) THEN     NbIMG2D=long(0) ELSE   NbIMG2D=long(NbIMG2D)

 PIOINFOIMG2DGRP=call_external(PIOLibIDLSO,'pioinfoimg2dgrp_tempoidl', $
        FLGname, $
        NbFLG, $
        IMG2Dtype, $
        IMG2Dname, $
        NbIMG2D, $
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
  IF (NbIMG2D GT 0) THEN BEGIN
     TP=IMG2Dtype
     TMPNAME=bytarr(128)
     IMG2Dtype=STRARR(NbIMG2D)
        FOR I=0L,NbIMG2D-1L DO BEGIN
            A=call_external(PIOLibIDLSO,'piogetstringidl',TP+I*128,TMPNAME)
            IMG2Dtype(I)=strmid(string(TMPNAME),0,A)
            END
    END
  IF (NbIMG2D GT 0) THEN BEGIN
     TP=IMG2Dname
     TMPNAME=bytarr(128)
     IMG2Dname=STRARR(NbIMG2D)
        FOR I=0L,NbIMG2D-1L DO BEGIN
            A=call_external(PIOLibIDLSO,'piogetstringidl',TP+I*128,TMPNAME)
            IMG2Dname(I)=strmid(string(TMPNAME),0,A)
            END
    END
    MESSAGE,'*',/CONTINUE
    MESSAGE,'NbIMG2D  :   '+STRING(NbIMG2D),/CONTINUE
    FOR I=0L,NbIMG2D-1L DO MESSAGE,'    '+'   '+STRING(IMG2DTYPE(I))+'   '+STRING(IMG2DNAME(I)) ,/CONTINUE
    MESSAGE,'*',/CONTINUE

 RETURN,PIOINFOIMG2DGRP
END
