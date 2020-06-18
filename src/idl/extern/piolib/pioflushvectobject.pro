
;--------------------------------------------------------------#
;
; AUTOMATICALLY GENERATED DO NOT MODIFY 

;--------------------------------------------------------------#




;	Wrapping PIOFlushVECTObject

FUNCTION PIOFLUSHVECTOBJECT,MyObject
     ON_ERROR,1
     PIOLibIDLSO=shared_lib_path('PIOLibIDL.so')

 PIOFLUSHVECTOBJECT=call_external(PIOLibIDLSO,'pioflushvectobject_tempoidl', $
        MyObject, $
               /L64_VALUE) 

 RETURN,PIOFLUSHVECTOBJECT
END
