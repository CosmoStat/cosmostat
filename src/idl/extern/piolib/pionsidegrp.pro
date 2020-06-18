
;--------------------------------------------------------------#
;
; AUTOMATICALLY GENERATED DO NOT MODIFY 

;--------------------------------------------------------------#




;	Wrapping PIONSideGrp

FUNCTION PIONSIDEGRP,MyGroup
     ON_ERROR,1
     PIOLibIDLSO=shared_lib_path('PIOLibIDL.so')

 PIONSIDEGRP=call_external(PIOLibIDLSO,'pionsidegrp_tempoidl', $
        MyGroup, $
               /L64_VALUE) 

 RETURN,PIONSIDEGRP
END
