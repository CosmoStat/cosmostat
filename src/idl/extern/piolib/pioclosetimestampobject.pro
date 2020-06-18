
;--------------------------------------------------------------#
;
; AUTOMATICALLY GENERATED DO NOT MODIFY 

;--------------------------------------------------------------#




;	Wrapping PIOCloseTIMESTAMPObject

FUNCTION PIOCLOSETIMESTAMPOBJECT,MyObject
     ON_ERROR,1
     PIOLibIDLSO=shared_lib_path('PIOLibIDL.so')
    MyObject=long64(MyObject)

 PIOCLOSETIMESTAMPOBJECT=call_external(PIOLibIDLSO,'pioclosetimestampobject_tempoidl', $
        MyObject, $
               /L64_VALUE) 

 RETURN,PIOCLOSETIMESTAMPOBJECT
END
