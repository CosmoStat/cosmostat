
;--------------------------------------------------------------#
;
; AUTOMATICALLY GENERATED DO NOT MODIFY 

;--------------------------------------------------------------#




;	Wrapping PIOCloseParFile

FUNCTION PIOCLOSEPARFILE,MyPar
     ON_ERROR,1
     PIOLibIDLSO=shared_lib_path('PIOLibIDL.so')
    MyPar=long64(MyPar)

 PIOCLOSEPARFILE=call_external(PIOLibIDLSO,'piocloseparfile_tempoidl', $
        MyPar, $
               /L64_VALUE) 

 RETURN,PIOCLOSEPARFILE
END
