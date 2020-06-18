FUNCTION PIOSETTABLE,TABLE,NAME,IMOUNIT

PIOLibIDLSO=shared_lib_path('PIOLibIDL.so')

tableunit=LONG64(0)

PIOERRMESS,call_external(PIOLibIDLSO,'pioinittable_tempoidl', $
              LONG64(TABLE.NBPOINT), PIOSTRING2C(TABLE.UNIT), $
              PIOSTRING2C(TABLE.NAME), PIOSTRING2C(TABLE.COMMENT), $
              PIOSTRING2C(TABLE.TYPE), TABLEUNIT,/L64_VALUE)

FOR I=0,TABLE.NBPOINT-1 DO BEGIN
    PIOERRMESS,call_external(PIOLibIDLSO,'piofilltable_tempoidl', $
              TABLE.X(I),TABLE.DX(I), $
              LONG64(I),TABLEUNIT,/L64_VALUE)
END

ERR = call_external(PIOLibIDLSO,'piosettableidl', $
              TABLEUNIT,PIOSTRING2C(NAME), $
              TABLE.TYPE,LONG64(IMOUNIT),/L64_VALUE)
PIOERRMESS,ERR
; Free memory allocated for table

return,call_external(PIOLibIDLSO,'piodeletetable_tempoidl',TABLEUNIT)

END
