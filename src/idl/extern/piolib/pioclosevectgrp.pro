
;--------------------------------------------------------------#
;
; AUTOMATICALLY GENERATED DO NOT MODIFY 

;--------------------------------------------------------------#




;	Wrapping PIOCloseVECTGrp

FUNCTION PIOCLOSEVECTGRP,MyGroup
     ON_ERROR,1
     PIOLibIDLSO=shared_lib_path('PIOLibIDL.so')
    MyGroup=long64(MyGroup)

 PIOCLOSEVECTGRP=call_external(PIOLibIDLSO,'pioclosevectgrp_tempoidl', $
        MyGroup, $
               /L64_VALUE) 

 RETURN,PIOCLOSEVECTGRP
END
