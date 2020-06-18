
;--------------------------------------------------------------#
;
; AUTOMATICALLY GENERATED DO NOT MODIFY 

;--------------------------------------------------------------#




;	Wrapping PIOSetEquinoxGrp

FUNCTION PIOSETEQUINOXGRP,Equinox,MyGroup
     ON_ERROR,1
     PIOLibIDLSO=shared_lib_path('PIOLibIDL.so')
    if (N_ELEMENTS(Equinox) EQ 0) THEN     Equinox=FLOAT(0) ELSE    Equinox=FLOAT(Equinox)
 
 PIOSETEQUINOXGRP=call_external(PIOLibIDLSO,'piosetequinoxgrp_tempoidl', $
        Equinox, $
        MyGroup, $
               /L64_VALUE) 

 RETURN,PIOSETEQUINOXGRP
END
