
;--------------------------------------------------------------#
;
; AUTOMATICALLY GENERATED DO NOT MODIFY 

;--------------------------------------------------------------#




;	Wrapping PIOGetCRValObject

FUNCTION PIOGETCRVALOBJECT,CRValX,CRValY,Objectname,MyGroup
     ON_ERROR,1
     PIOLibIDLSO=shared_lib_path('PIOLibIDL.so')
    if (N_ELEMENTS(CRValX) EQ 0) THEN     CRValX=DOUBLE(0) ELSE    CRValX=DOUBLE(CRValX)
     if (N_ELEMENTS(CRValY) EQ 0) THEN     CRValY=DOUBLE(0) ELSE    CRValY=DOUBLE(CRValY)
     Objectname_TMP=BYTARR(128)
    Objectname_TMP(*)=0
    if (N_ELEMENTS(Objectname) GT 0) THEN if (STRLEN(Objectname) GT 0) THEN Objectname_TMP(0:STRLEN(Objectname)-1)=BYTE(Objectname)
IF (N_PARAMS() LT 4) THEN BEGIN
    MyGroup=long64(0)
END ELSE BEGIN
    MyGroup=long64(MyGroup)
END

 PIOGETCRVALOBJECT=call_external(PIOLibIDLSO,'piogetcrvalobject_tempoidl', $
        CRValX, $
        CRValY, $
        Objectname_TMP, $
        MyGroup, $
               /L64_VALUE) 
  Objectname=STRING(Objectname_TMP)

 RETURN,PIOGETCRVALOBJECT
END
