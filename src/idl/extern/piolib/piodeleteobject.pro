
;--------------------------------------------------------------#
;
; AUTOMATICALLY GENERATED DO NOT MODIFY 

;--------------------------------------------------------------#




;	Wrapping PIODeleteObject

FUNCTION PIODELETEOBJECT,Objectname,MyGroup
     ON_ERROR,1
     PIOLibIDLSO=shared_lib_path('PIOLibIDL.so')
    Objectname_TMP=BYTARR(128)
    Objectname_TMP(*)=0
    if (N_ELEMENTS(Objectname) GT 0) THEN if (STRLEN(Objectname) GT 0) THEN Objectname_TMP(0:STRLEN(Objectname)-1)=BYTE(Objectname)
IF (N_PARAMS() LT 2) THEN BEGIN
    MyGroup=long64(0)
END ELSE BEGIN
    MyGroup=long64(MyGroup)
END

 PIODELETEOBJECT=call_external(PIOLibIDLSO,'piodeleteobject_tempoidl', $
        Objectname_TMP, $
        MyGroup, $
               /L64_VALUE) 

 RETURN,PIODELETEOBJECT
END
