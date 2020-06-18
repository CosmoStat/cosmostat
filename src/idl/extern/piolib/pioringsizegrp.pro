
;--------------------------------------------------------------#
;
; AUTOMATICALLY GENERATED DO NOT MODIFY 

;--------------------------------------------------------------#




;	Wrapping PIORingSizeGrp

FUNCTION PIORINGSIZEGRP,MyGroup
     ON_ERROR,1
     PIOLibIDLSO=shared_lib_path('PIOLibIDL.so')

 PIORINGSIZEGRP=call_external(PIOLibIDLSO,'pioringsizegrp_tempoidl', $
        MyGroup, $
               /L64_VALUE) 

 RETURN,PIORINGSIZEGRP
END
