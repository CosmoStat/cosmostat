
;--------------------------------------------------------------#
;
; AUTOMATICALLY GENERATED DO NOT MODIFY 

;--------------------------------------------------------------#




;	Wrapping PIODeleteTable2col

FUNCTION PIODELETETABLE2COL,MyTable2col
     ON_ERROR,1
     PIOLibIDLSO=shared_lib_path('PIOLibIDL.so')
 IF (N_ELEMENTS(MyTable2col) EQ 0) THEN MyTable2col=LONG64(0) ELSE MyTable2col=LONG64(MyTable2col)

 PIODELETETABLE2COL=call_external(PIOLibIDLSO,'piodeletetable2col_tempoidl', $
        MyTable2col, $
               /L64_VALUE) 

 RETURN,PIODELETETABLE2COL
END
