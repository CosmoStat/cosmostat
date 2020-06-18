
;--------------------------------------------------------------#
;
; AUTOMATICALLY GENERATED DO NOT MODIFY 

;--------------------------------------------------------------#




;	Wrapping PIOCreateMAPGrp

FUNCTION PIOCREATEMAPGRP,Groupname,Coordsys,Ordering,NSide
     ON_ERROR,1
     PIOLibIDLSO=shared_lib_path('PIOLibIDL.so')
    Groupname_TMP=BYTARR(128)
    Groupname_TMP(*)=0
    if (N_ELEMENTS(Groupname) GT 0) THEN if (STRLEN(Groupname) GT 0) THEN Groupname_TMP(0:STRLEN(Groupname)-1)=BYTE(Groupname)
    Coordsys_TMP=BYTARR(128)
    Coordsys_TMP(*)=0
    if (N_ELEMENTS(Coordsys) GT 0) THEN if (STRLEN(Coordsys) GT 0) THEN Coordsys_TMP(0:STRLEN(Coordsys)-1)=BYTE(Coordsys)
    Ordering_TMP=BYTARR(128)
    Ordering_TMP(*)=0
    if (N_ELEMENTS(Ordering) GT 0) THEN if (STRLEN(Ordering) GT 0) THEN Ordering_TMP(0:STRLEN(Ordering)-1)=BYTE(Ordering)
    if (N_ELEMENTS(NSide) EQ 0) THEN     NSide=LONG64(0) ELSE    NSide=LONG64(NSide)
 
 PIOCREATEMAPGRP=call_external(PIOLibIDLSO,'piocreatemapgrp_tempoidl', $
        Groupname_TMP, $
        Coordsys_TMP, $
        Ordering_TMP, $
        NSide, $
               /L64_VALUE) 

 RETURN,PIOCREATEMAPGRP
END
