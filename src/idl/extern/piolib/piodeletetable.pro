
;--------------------------------------------------------------#
;
; AUTOMATICALLY GENERATED DO NOT MODIFY 

;--------------------------------------------------------------#




;	Wrapping PIODeleteTable

FUNCTION PIODELETETABLE,MyTable
     ON_ERROR,1
     PIOLibIDLSO=shared_lib_path('PIOLibIDL.so')
 IF (N_ELEMENTS(MyTable) EQ 0) THEN MyTable=LONG64(0) ELSE MyTable=LONG64(MyTable)

 PIODELETETABLE=call_external(PIOLibIDLSO,'piodeletetable_tempoidl', $
        MyTable, $
               /L64_VALUE) 

 RETURN,PIODELETETABLE
END
