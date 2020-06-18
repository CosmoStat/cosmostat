
;--------------------------------------------------------------#
;
; AUTOMATICALLY GENERATED DO NOT MODIFY 

;--------------------------------------------------------------#




;	Wrapping PIOCreateTOIObject

FUNCTION PIOCREATETOIOBJECT,Objectname,typein,MyGroup
     ON_ERROR,1
     PIOLibIDLSO=shared_lib_path('PIOLibIDL.so')
    Objectname_TMP=BYTARR(128)
    Objectname_TMP(*)=0
    if (N_ELEMENTS(Objectname) GT 0) THEN if (STRLEN(Objectname) GT 0) THEN Objectname_TMP(0:STRLEN(Objectname)-1)=BYTE(Objectname)
    typein_TMP=BYTARR(128)
    typein_TMP(*)=0
    if (N_ELEMENTS(typein) GT 0) THEN if (STRLEN(typein) GT 0) THEN typein_TMP(0:STRLEN(typein)-1)=BYTE(typein)
IF (N_PARAMS() LT 3) THEN BEGIN
    MyGroup=long64(0)
END ELSE BEGIN
    MyGroup=long64(MyGroup)
END

 PIOCREATETOIOBJECT=call_external(PIOLibIDLSO,'piocreatetoiobject_tempoidl', $
        Objectname_TMP, $
        typein_TMP, $
        MyGroup, $
               /L64_VALUE) 

 RETURN,PIOCREATETOIOBJECT
END
