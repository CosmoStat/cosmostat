
;--------------------------------------------------------------#
;
; AUTOMATICALLY GENERATED DO NOT MODIFY 

;--------------------------------------------------------------#




;	Wrapping PIOGetDatabaseName

FUNCTION PIOGETDATABASENAME,Database,HTMLPath
     ON_ERROR,1
     PIOLibIDLSO=shared_lib_path('PIOLibIDL.so')
    Database_TMP=BYTARR(128)
    Database_TMP(*)=0
    if (N_ELEMENTS(Database) GT 0) THEN if (STRLEN(Database) GT 0) THEN Database_TMP(0:STRLEN(Database)-1)=BYTE(Database)
    HTMLPath_TMP=BYTARR(128)
    HTMLPath_TMP(*)=0
    if (N_ELEMENTS(HTMLPath) GT 0) THEN if (STRLEN(HTMLPath) GT 0) THEN HTMLPath_TMP(0:STRLEN(HTMLPath)-1)=BYTE(HTMLPath)

 PIOGETDATABASENAME=call_external(PIOLibIDLSO,'piogetdatabasename_tempoidl', $
        Database_TMP, $
        HTMLPath_TMP, $
               /L64_VALUE) 
  Database=STRING(Database_TMP)
  HTMLPath=STRING(HTMLPath_TMP)

 RETURN,PIOGETDATABASENAME
END
