
;--------------------------------------------------------------#
;
; AUTOMATICALLY GENERATED DO NOT MODIFY 

;--------------------------------------------------------------#




;	Wrapping PIOGetVersion

FUNCTION PIOGETVERSION,VersionTag
     ON_ERROR,1
     PIOLibIDLSO=shared_lib_path('PIOLibIDL.so')
    VersionTag_TMP=BYTARR(128)
    VersionTag_TMP(*)=0
    if (N_ELEMENTS(VersionTag) GT 0) THEN if (STRLEN(VersionTag) GT 0) THEN VersionTag_TMP(0:STRLEN(VersionTag)-1)=BYTE(VersionTag)

 PIOGETVERSION=call_external(PIOLibIDLSO,'piogetversion_tempoidl', $
        VersionTag_TMP, $
               /L64_VALUE) 
  VersionTag=STRING(VersionTag_TMP)

 RETURN,PIOGETVERSION
END
