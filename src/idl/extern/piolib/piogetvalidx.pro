
;--------------------------------------------------------------#
;
; AUTOMATICALLY GENERATED DO NOT MODIFY 

;--------------------------------------------------------------#




;	Wrapping PIOGetValIdx

FUNCTION PIOGETVALIDX,Value,Obj,MyGroup
     ON_ERROR,1
     PIOLibIDLSO=shared_lib_path('PIOLibIDL.so')
    if (N_ELEMENTS(Value) EQ 0) THEN     Value=DOUBLE(0) ELSE    Value=DOUBLE(Value)
     Obj_TMP=BYTARR(128)
    Obj_TMP(*)=0
    if (N_ELEMENTS(Obj) GT 0) THEN if (STRLEN(Obj) GT 0) THEN Obj_TMP(0:STRLEN(Obj)-1)=BYTE(Obj)
IF (N_PARAMS() LT 3) THEN BEGIN
    MyGroup=long64(0)
END ELSE BEGIN
    MyGroup=long64(MyGroup)
END

 PIOGETVALIDX=call_external(PIOLibIDLSO,'piogetvalidx_tempoidl', $
        Value, $
        Obj_TMP, $
        MyGroup, $
               /L64_VALUE) 

 RETURN,PIOGETVALIDX
END
