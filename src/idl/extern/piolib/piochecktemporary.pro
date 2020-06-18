
;--------------------------------------------------------------#
;
; AUTOMATICALLY GENERATED DO NOT MODIFY 

;--------------------------------------------------------------#




;	Wrapping PIOCheckTemporary

FUNCTION PIOCHECKTEMPORARY,ObjectName
     ON_ERROR,1
     PIOLibIDLSO=shared_lib_path('PIOLibIDL.so')
    ObjectName_TMP=BYTARR(128)
    ObjectName_TMP(*)=0
    if (N_ELEMENTS(ObjectName) GT 0) THEN if (STRLEN(ObjectName) GT 0) THEN ObjectName_TMP(0:STRLEN(ObjectName)-1)=BYTE(ObjectName)

 PIOCHECKTEMPORARY=call_external(PIOLibIDLSO,'piochecktemporary_tempoidl', $
        ObjectName_TMP, $
               /L64_VALUE) 

 RETURN,PIOCHECKTEMPORARY
END
