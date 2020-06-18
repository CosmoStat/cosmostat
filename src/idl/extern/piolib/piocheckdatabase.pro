
;--------------------------------------------------------------#
;
; AUTOMATICALLY GENERATED DO NOT MODIFY 

;--------------------------------------------------------------#




;	Wrapping PIOCheckDataBase

FUNCTION PIOCHECKDATABASE,Database
     ON_ERROR,1
     PIOLibIDLSO=shared_lib_path('PIOLibIDL.so')
    Database_TMP=BYTARR(128)
    Database_TMP(*)=0
    if (N_ELEMENTS(Database) GT 0) THEN if (STRLEN(Database) GT 0) THEN Database_TMP(0:STRLEN(Database)-1)=BYTE(Database)

 PIOCHECKDATABASE=call_external(PIOLibIDLSO,'piocheckdatabase_tempoidl', $
        Database_TMP, $
               /L64_VALUE) 

 RETURN,PIOCHECKDATABASE
END
