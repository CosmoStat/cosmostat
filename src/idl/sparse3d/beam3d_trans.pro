;+
; NAME: 
;       BEAM3D_TRANS
;
; PURPOSE:
;       Compute the 3D beamlet  transform of a cube. The output is an IDL structure.
;       A band at scale j (j=0..NBRSCALE-1) can be extracted using the function
;       function rid3d_getband(Rid, j) (ex: Scale2 = rid3d_getband(BeamTrans, 2))
;       and a band can be inserted in the transformation using the routine  rid3d_putband
;       (ex:  rid3d_putband, BeamTrans, Scale2, 2).
;
; CALLING:
;     rid3d_trans, Cube, Trans, Opt=Opt 
;
; INPUTS:
;     Cube -- IDL 3D array: Input cube be transformed 
;     
; OUTPUTS:
;     Trans -- IDL structures with the following fields:  
;         NBRSCALE  -- LONG: Number of scales of the ridgelet transform
;         NXBeam     -- LONG: Number of pixels in the x-axis direction
;         NYBeam     -- LONG: Number of pixels in the y-axis direction
;         NZBeam     -- LONG: Number of pixels in the z-axis direction
;         COEF      -- 3D IDL array: Beamlet coefficients
;                             cube which contains all ridgelet coefficients
;         BSIZE     --  LONG: Block size used in the ridgelet transform
;         NXB       -- LONG: Number of blocks in the x-axis direction
;         NYB       -- LONG: Number of blocks in the y-axis direction
;         NZB       -- LONG: Number of blocks in the z-axis direction
;         OVERLAP   -- LONG: is equal to 1 if blocks are overlapping
;         HEADTRANS -- IDL STRING Array: contains the Fits Header of the decomposition
;
; KEYWORDS:
;      OPT: string which contains the differents options. Options are:
;         [-b BlockSize]
;             Block Size. Default is image size. 
;         [-O]
;             Block overlapping. Default is no. 
;         [-v]
;             Verbose. Default is no.
;
; EXTERNAL CALLS
;           beam3d_trans (C++ program)
;
; EXAMPLE:
;     beam3d_trans, Cube, Rid
;       beamlet transform with default options.
;
; HISTORY:
;       Written: Jean-Luc Starck 2005.
;       September, 2005 File creation
;-
;-----------------------------------------------------------------
 
pro beam3d_trans, Ima, BeamTrans, OPT=OPT 

if N_PARAMS() LT 2 then begin 
        spawn, 'beam3d_trans'
        print, 'CALL SEQUENCE: beam3d_trans, Signal, Struct_Out, OPT=Opt'
        goto, DONE
        end

Nx = (size(Ima))[1]
Ny = (size(Ima))[2]

NameIma = 'xx_signal.fits'
NameResult = 'xx_result'
NameResultFits = 'xx_result.bet'
; if TypeOrder  eq 'RING    '  or  TypeOrder  eq 'ring    ' then Order=1

writefits,  NameIma,  Ima
if not keyword_set(OPT) then OPT = " "
 
com = 'beam3d_trans  ' + OPT + ' '+ NameIma  + ' ' +  NameResult
spawn, com

Beam = readfits(NameResultFits, HeadTrans, /silent); 
;help, Beam
;help, HeadTrans
NxBeam = (size(Beam))[1]
NyBeam = (size(Beam))[2]
NzBeam = (size(Beam))[3]
Nx = FXPAR(HeadTrans, "NX")
NY = FXPAR(HeadTrans, "NY")
NZ = FXPAR(HeadTrans, "NZ")

Bsize = FXPAR( HeadTrans, "BLOCKSIZ")
Nxb = FXPAR(HeadTrans, "NXB")
Nyb = FXPAR(HeadTrans, "NYB")
Nzb = FXPAR(HeadTrans, "NZB")

Overlap = FXPAR(HeadTrans, "OVERLAP")
NbrScale = FXPAR(HeadTrans, "NBRSCALE")
TypeTrans = FXPAR(HeadTrans, "TYPE_TRA")                               

;print, ' Nl = ' , Nl, ' Nc = ', Nc, ', Bsize =', Bsize, ' NbrScale =  ' , NbrScale
;print, ' Nlb = ' , Nlb, ' Ncb = ', Ncb, ' Overlap = ', Overlap, ' TypeTrans = ', TypeTrans
;info, Rid-ImaRid
;  mr1d_struc, result, info, output, band=band

; TabNorm=[0.8546,0.215, 0.1265,0.097,0.0669,0.0483,0.0342,0.0242,0.0171,0.0121]
BeamTrans = {NbrScale : NbrScale, NxBeam : NxBeam, NyBeam:NyBeam,  NzBeam : NzBeam, Coef : Beam, Bsize : Bsize, Nxb: Nxb, Nyb:Nyb, Nzb: Nzb, $
            Overlap: Overlap, HeadTrans:HeadTrans}
DONE:

end

 

pro trid3d, e,r, opt=opt
N=8
e = fltarr(N,N,N)
e(*,N/2,*)=1
e = readfits('b.fits')
rid3d_trans, e, BeamTrans, OPT=OPT 
 
NBeamTrans = BeamTrans
NBeamTrans.coef[*] = 0
for j=0,BeamTrans.NbrScale-1 do begin
   B = rid3d_getband(BeamTrans, j)
   rid3d_putband,NBeamTrans, B, j
   end

info, NBeamTrans.coef - BeamTrans.coef 
rid3d_rec, BeamTrans, r
end
