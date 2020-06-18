
;--------------------------------------------------------------#
;
; AUTOMATICALLY GENERATED DO NOT MODIFY 

;--------------------------------------------------------------#




;	Wrapping PIOCreateIMG3DGrp

FUNCTION PIOCREATEIMG3DGRP,GroupName,Proj,Coordsys,NAXIS1,NAXIS2,NAXIS3,PixSize1,PixSize2,CRPixX,CRPixY,PV
     ON_ERROR,1
     PIOLibIDLSO=shared_lib_path('PIOLibIDL.so')
    GroupName_TMP=BYTARR(128)
    GroupName_TMP(*)=0
    if (N_ELEMENTS(GroupName) GT 0) THEN if (STRLEN(GroupName) GT 0) THEN GroupName_TMP(0:STRLEN(GroupName)-1)=BYTE(GroupName)
    Proj_TMP=BYTARR(128)
    Proj_TMP(*)=0
    if (N_ELEMENTS(Proj) GT 0) THEN if (STRLEN(Proj) GT 0) THEN Proj_TMP(0:STRLEN(Proj)-1)=BYTE(Proj)
    Coordsys_TMP=BYTARR(128)
    Coordsys_TMP(*)=0
    if (N_ELEMENTS(Coordsys) GT 0) THEN if (STRLEN(Coordsys) GT 0) THEN Coordsys_TMP(0:STRLEN(Coordsys)-1)=BYTE(Coordsys)
    if (N_ELEMENTS(NAXIS1) EQ 0) THEN     NAXIS1=LONG64(0) ELSE    NAXIS1=LONG64(NAXIS1)
     if (N_ELEMENTS(NAXIS2) EQ 0) THEN     NAXIS2=LONG64(0) ELSE    NAXIS2=LONG64(NAXIS2)
     if (N_ELEMENTS(NAXIS3) EQ 0) THEN     NAXIS3=LONG64(0) ELSE    NAXIS3=LONG64(NAXIS3)
     if (N_ELEMENTS(PixSize1) EQ 0) THEN     PixSize1=FLOAT(0) ELSE    PixSize1=FLOAT(PixSize1)
     if (N_ELEMENTS(PixSize2) EQ 0) THEN     PixSize2=FLOAT(0) ELSE    PixSize2=FLOAT(PixSize2)
     if (N_ELEMENTS(CRPixX) EQ 0) THEN     CRPixX=DOUBLE(0) ELSE    CRPixX=DOUBLE(CRPixX)
     if (N_ELEMENTS(CRPixY) EQ 0) THEN     CRPixY=DOUBLE(0) ELSE    CRPixY=DOUBLE(CRPixY)
 IF (N_PARAMS() LT 11) THEN BEGIN
    PV=long64(0)
END ELSE BEGIN
    PV=long64(PV)
END

 PIOCREATEIMG3DGRP=call_external(PIOLibIDLSO,'piocreateimg3dgrp_tempoidl', $
        GroupName_TMP, $
        Proj_TMP, $
        Coordsys_TMP, $
        NAXIS1, $
        NAXIS2, $
        NAXIS3, $
        PixSize1, $
        PixSize2, $
        CRPixX, $
        CRPixY, $
        PV, $
               /L64_VALUE) 

 RETURN,PIOCREATEIMG3DGRP
END
