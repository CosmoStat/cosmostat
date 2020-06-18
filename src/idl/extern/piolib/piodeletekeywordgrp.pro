
;--------------------------------------------------------------#
;
; AUTOMATICALLY GENERATED DO NOT MODIFY 

;--------------------------------------------------------------#




;	Wrapping PIODeleteKeywordGrp

FUNCTION PIODELETEKEYWORDGRP,InKeyword,MyGroup
     ON_ERROR,1
     PIOLibIDLSO=shared_lib_path('PIOLibIDL.so')
    InKeyword_TMP=BYTARR(128)
    InKeyword_TMP(*)=0
    if (N_ELEMENTS(InKeyword) GT 0) THEN if (STRLEN(InKeyword) GT 0) THEN InKeyword_TMP(0:STRLEN(InKeyword)-1)=BYTE(InKeyword)

 PIODELETEKEYWORDGRP=call_external(PIOLibIDLSO,'piodeletekeywordgrp_tempoidl', $
        InKeyword_TMP, $
        MyGroup, $
               /L64_VALUE) 

 RETURN,PIODELETEKEYWORDGRP
END
