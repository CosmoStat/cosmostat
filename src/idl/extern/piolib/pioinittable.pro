FUNCTION PIOINITTABLE,NBPOINT,UNIT,NAME,COMMENT,TYPE

PIOLibIDLSO=shared_lib_path('PIOLibIDL.so')

TABLEUNIT=LONG64(0)

PIOERRMESS,call_external(PIOLibIDLSO,'pioinittable_tempoidl', $
              LONG64(NBPOINT),PIOSTRING2C(UNIT),PIOSTRING2C(NAME), $
              PIOSTRING2C(COMMENT),PIOSTRING2C(TYPE),TABLEUNIT,/L64_VALUE)

CASE (STRING(TYPE)) OF
        "PIOLONG": BEGIN
            X=LONG64(0)
            DATAX = LON64ARR(NBPOINT)
            END
        "PIOBYTE": BEGIN
            X=BYTE(0)
            DATAX = BYTARR(NBPOINT)
            END
        "PIOSHORT": BEGIN
            X=0
            DATAX = INTARR(NBPOINT)
            END
        "PIOFLOAT": BEGIN
            X=FLOAT(0)
            DATAX = FLTARR(NBPOINT)
            END
        "PIODOUBLE": BEGIN
            X=DOUBLE(0)
            DATAX = DBLARR(NBPOINT)
            END
        "PIOCOMPLEX": BEGIN
            X=COMPLEX(0,0)
            DATAX = COMPLEXARR(NBPOINT)
            END
        "PIODBLCOMPLEX": BEGIN
            X=DCOMPLEX(0,0)
            DATAX = DCOMPLEXARR(NBPOINT)
            END
        ELSE: return,-1
    END

TABLE ={NBPOINT:NBPOINT, $
        UNIT:UNIT,NAME:NAME, $
        COMMENT:COMMENT, TYPE:TYPE, $
        TABLEUNIT:TABLEUNIT, $
        X:DATAX,DX:DATAX}

; Free memory allocated for table

PIOERRMESS,call_external(PIOLibIDLSO,'piodeletetable_tempoidl', $
                        LONG64(TABLEUNIT))

RETURN,TABLE
END
