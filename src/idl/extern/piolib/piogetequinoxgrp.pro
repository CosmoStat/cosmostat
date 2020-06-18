
;--------------------------------------------------------------#
;
; AUTOMATICALLY GENERATED DO NOT MODIFY 

;--------------------------------------------------------------#




;	Wrapping PIOGetEquinoxGrp

FUNCTION PIOGETEQUINOXGRP,Equinox,MyGroup
     ON_ERROR,1
     PIOLibIDLSO=shared_lib_path('PIOLibIDL.so')
    if (N_ELEMENTS(Equinox) EQ 0) THEN     Equinox=FLOAT(0) ELSE    Equinox=FLOAT(Equinox)
 
 PIOGETEQUINOXGRP=call_external(PIOLibIDLSO,'piogetequinoxgrp_tempoidl', $
        Equinox, $
        MyGroup, $
               /L64_VALUE) 

 RETURN,PIOGETEQUINOXGRP
END
