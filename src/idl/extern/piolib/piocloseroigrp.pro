
;--------------------------------------------------------------#
;
; AUTOMATICALLY GENERATED DO NOT MODIFY 

;--------------------------------------------------------------#




;	Wrapping PIOCloseROIGrp

FUNCTION PIOCLOSEROIGRP,object
     ON_ERROR,1
     PIOLibIDLSO=shared_lib_path('PIOLibIDL.so')
    object=long64(object)

 PIOCLOSEROIGRP=call_external(PIOLibIDLSO,'piocloseroigrp_tempoidl', $
        object, $
               /L64_VALUE) 

 RETURN,PIOCLOSEROIGRP
END
