
;--------------------------------------------------------------#
;
; AUTOMATICALLY GENERATED DO NOT MODIFY 

;--------------------------------------------------------------#




;	Wrapping PIOCreateIMG3DObject

FUNCTION PIOCREATEIMG3DOBJECT,Objectname,typein,CRValX,CRValY,MyGroup
     ON_ERROR,1
     PIOLibIDLSO=shared_lib_path('PIOLibIDL.so')
    Objectname_TMP=BYTARR(128)
    Objectname_TMP(*)=0
    if (N_ELEMENTS(Objectname) GT 0) THEN if (STRLEN(Objectname) GT 0) THEN Objectname_TMP(0:STRLEN(Objectname)-1)=BYTE(Objectname)
    typein_TMP=BYTARR(128)
    typein_TMP(*)=0
    if (N_ELEMENTS(typein) GT 0) THEN if (STRLEN(typein) GT 0) THEN typein_TMP(0:STRLEN(typein)-1)=BYTE(typein)
    if (N_ELEMENTS(CRValX) EQ 0) THEN     CRValX=DOUBLE(0) ELSE    CRValX=DOUBLE(CRValX)
     if (N_ELEMENTS(CRValY) EQ 0) THEN     CRValY=DOUBLE(0) ELSE    CRValY=DOUBLE(CRValY)
 IF (N_PARAMS() LT 5) THEN BEGIN
    MyGroup=long64(0)
END ELSE BEGIN
    MyGroup=long64(MyGroup)
END

 PIOCREATEIMG3DOBJECT=call_external(PIOLibIDLSO,'piocreateimg3dobject_tempoidl', $
        Objectname_TMP, $
        typein_TMP, $
        CRValX, $
        CRValY, $
        MyGroup, $
               /L64_VALUE) 

 RETURN,PIOCREATEIMG3DOBJECT
END
