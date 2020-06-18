
;--------------------------------------------------------------#
;
; AUTOMATICALLY GENERATED DO NOT MODIFY 

;--------------------------------------------------------------#




;	Wrapping PIOInfoIMG3DGrp

FUNCTION PIOINFOIMG3DGRP,MyGroup,FLGname=FLGname,NbFLG=NbFLG,IMG3Dtype=IMG3Dtype,IMG3Dname=IMG3Dname,NbIMG3D=NbIMG3D
     ON_ERROR,1
     PIOLibIDLSO=shared_lib_path('PIOLibIDL.so')
    FLGname=long64(0)
    if (N_ELEMENTS(NbFLG) EQ 0) THEN     NbFLG=long(0) ELSE   NbFLG=long(NbFLG)
    IMG3Dtype=long64(0)
    IMG3Dname=long64(0)
    if (N_ELEMENTS(NbIMG3D) EQ 0) THEN     NbIMG3D=long(0) ELSE   NbIMG3D=long(NbIMG3D)

 PIOINFOIMG3DGRP=call_external(PIOLibIDLSO,'pioinfoimg3dgrp_tempoidl', $
        FLGname, $
        NbFLG, $
        IMG3Dtype, $
        IMG3Dname, $
        NbIMG3D, $
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
  IF (NbIMG3D GT 0) THEN BEGIN
     TP=IMG3Dtype
     TMPNAME=bytarr(128)
     IMG3Dtype=STRARR(NbIMG3D)
        FOR I=0L,NbIMG3D-1L DO BEGIN
            A=call_external(PIOLibIDLSO,'piogetstringidl',TP+I*128,TMPNAME)
            IMG3Dtype(I)=strmid(string(TMPNAME),0,A)
            END
    END
  IF (NbIMG3D GT 0) THEN BEGIN
     TP=IMG3Dname
     TMPNAME=bytarr(128)
     IMG3Dname=STRARR(NbIMG3D)
        FOR I=0L,NbIMG3D-1L DO BEGIN
            A=call_external(PIOLibIDLSO,'piogetstringidl',TP+I*128,TMPNAME)
            IMG3Dname(I)=strmid(string(TMPNAME),0,A)
            END
    END
    MESSAGE,'*',/CONTINUE
    MESSAGE,'NbIMG3D  :   '+STRING(NbIMG3D),/CONTINUE
    FOR I=0L,NbIMG3D-1L DO MESSAGE,'    '+'   '+STRING(IMG3DTYPE(I))+'   '+STRING(IMG3DNAME(I)) ,/CONTINUE
    MESSAGE,'*',/CONTINUE

 RETURN,PIOINFOIMG3DGRP
END
