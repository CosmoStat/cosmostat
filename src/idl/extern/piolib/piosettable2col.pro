FUNCTION PIOSETTABLE2COL,TABLE,NAME,IMOUNIT

PIOLibIDLSO=shared_lib_path('PIOLibIDL.so')


tableunit=LONG64(0)
PIOERRMESS,call_external(PIOLibIDLSO,'pioinittable2col_tempoidl', $
              LONG64(TABLE.NBPOINT), $
              PIOSTRING2C(TABLE.XUNIT),PIOSTRING2C(TABLE.YUNIT), $
              PIOSTRING2C(TABLE.XNAME),PIOSTRING2C(TABLE.YNAME), $
                         PIOSTRING2C(TABLE.COMMENT), $
              PIOSTRING2C(TABLE.XTYPE),PIOSTRING2C(TABLE.YTYPE), $
                         TABLEUNIT,/L64_VALUE)

FOR I=0,TABLE.NBPOINT-1 DO BEGIN
    PIOERRMESS,call_external(PIOLibIDLSO,'piofilltable2col_tempoidl', $
              TABLE.X(I),TABLE.DX(I),TABLE.Y(i),TABLE.DY(I), $
              LONG64(I),TABLEUNIT,/L64_VALUE)
END

ERR = call_external(PIOLibIDLSO,'piosettable2col_tempoidl', $
              TABLEUNIT,PIOSTRING2C(NAME), $
              TABLE.XTYPE,TABLE.YTYPE,LONG64(IMOUNIT),/L64_VALUE)
PIOERRMESS,ERR
; Free memory allocated for table

return,call_external(PIOLibIDLSO,'piodeletetable2col_tempoidl',TABLEUNIT)

END
