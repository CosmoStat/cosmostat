
;--------------------------------------------------------------#
;
; AUTOMATICALLY GENERATED DO NOT MODIFY 

;--------------------------------------------------------------#




;	Wrapping PIOGetEndObjectIdx

FUNCTION PIOGETENDOBJECTIDX,ObjectName,MyGroup
     ON_ERROR,1
     PIOLibIDLSO=shared_lib_path('PIOLibIDL.so')
    ObjectName_TMP=BYTARR(128)
    ObjectName_TMP(*)=0
    if (N_ELEMENTS(ObjectName) GT 0) THEN if (STRLEN(ObjectName) GT 0) THEN ObjectName_TMP(0:STRLEN(ObjectName)-1)=BYTE(ObjectName)
IF (N_PARAMS() LT 2) THEN BEGIN
    MyGroup=long64(0)
END ELSE BEGIN
    MyGroup=long64(MyGroup)
END

 PIOGETENDOBJECTIDX=call_external(PIOLibIDLSO,'piogetendobjectidx_tempoidl', $
        ObjectName_TMP, $
        MyGroup, $
               /L64_VALUE) 

 RETURN,PIOGETENDOBJECTIDX
END
