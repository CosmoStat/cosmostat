
;--------------------------------------------------------------#
;
; AUTOMATICALLY GENERATED DO NOT MODIFY 

;--------------------------------------------------------------#




;	Wrapping PIOGetTAB3DDeep

FUNCTION PIOGETTAB3DDEEP,FirstAx,LastAx,ObjectName,MyGroup
     ON_ERROR,1
     PIOLibIDLSO=shared_lib_path('PIOLibIDL.so')
    if (N_ELEMENTS(FirstAx) EQ 0) THEN     FirstAx=LONG64(0) ELSE    FirstAx=LONG64(FirstAx)
     if (N_ELEMENTS(LastAx) EQ 0) THEN     LastAx=LONG64(0) ELSE    LastAx=LONG64(LastAx)
     ObjectName_TMP=BYTARR(128)
    ObjectName_TMP(*)=0
    if (N_ELEMENTS(ObjectName) GT 0) THEN if (STRLEN(ObjectName) GT 0) THEN ObjectName_TMP(0:STRLEN(ObjectName)-1)=BYTE(ObjectName)

 PIOGETTAB3DDEEP=call_external(PIOLibIDLSO,'piogettab3ddeep_tempoidl', $
        FirstAx, $
        LastAx, $
        ObjectName_TMP, $
        MyGroup, $
               /L64_VALUE) 
  ObjectName=STRING(ObjectName_TMP)

 RETURN,PIOGETTAB3DDEEP
END
