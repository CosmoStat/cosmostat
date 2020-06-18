
;--------------------------------------------------------------#
;
; AUTOMATICALLY GENERATED DO NOT MODIFY 

;--------------------------------------------------------------#




;	Wrapping PIOInfoIMO

FUNCTION PIOINFOIMO,MajorVersion,MinorVersion,MyGroup
     ON_ERROR,1
     PIOLibIDLSO=shared_lib_path('PIOLibIDL.so')
    if (N_ELEMENTS(MajorVersion) EQ 0) THEN     MajorVersion=LONG64(0) ELSE    MajorVersion=LONG64(MajorVersion)
     if (N_ELEMENTS(MinorVersion) EQ 0) THEN     MinorVersion=LONG64(0) ELSE    MinorVersion=LONG64(MinorVersion)
 
 PIOINFOIMO=call_external(PIOLibIDLSO,'pioinfoimo_tempoidl', $
        MajorVersion, $
        MinorVersion, $
        MyGroup, $
               /L64_VALUE) 

 RETURN,PIOINFOIMO
END
