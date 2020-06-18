
;--------------------------------------------------------------#
;
; AUTOMATICALLY GENERATED DO NOT MODIFY 

;--------------------------------------------------------------#




;	Wrapping PIOGetLatpoleGrp

FUNCTION PIOGETLATPOLEGRP,LatPole,MyGroup
     ON_ERROR,1
     PIOLibIDLSO=shared_lib_path('PIOLibIDL.so')
    if (N_ELEMENTS(LatPole) EQ 0) THEN     LatPole=DOUBLE(0) ELSE    LatPole=DOUBLE(LatPole)
 
 PIOGETLATPOLEGRP=call_external(PIOLibIDLSO,'piogetlatpolegrp_tempoidl', $
        LatPole, $
        MyGroup, $
               /L64_VALUE) 

 RETURN,PIOGETLATPOLEGRP
END
