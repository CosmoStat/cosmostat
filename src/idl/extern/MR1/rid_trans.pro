;+
; NAME: 
;       RID_TRANS
;
; PURPOSE:
;       Compute the ridgelet  transform of an image. The output is a IDL structure.
;       A band at scale j (j=0..NBRSCALE-1) can be extracted using the function
;       function rid_getband(Rid, j) (ex: Scale2 = rid_getband(RidTrans, 2))
;       and a band can be inserted in the transformation using the routine  rid_putband
;       (ex:  rid_putband, RidTrans, Scale2, 2).
;
; CALLING:
;     rid_trans, Imag, Trans, Opt=Opt 
;
; INPUTS:
;     Imag -- IDL 2D array: Input image be transformed 
;     
; OUTPUTS:
;     Trans -- IDL structures with the following fields:  
;         NBRSCALE  -- LONG: Number of scales of the ridgelet transform
;         NXRID     -- LONG: Number of pixels in the x-axis direction
;         NYRID     -- LONG: Number of pixels in the y-axis direction
;         COEF      -- 2D IDL array: Ridgelet coefficients
;                             Image which contains all ridgelet coefficients
;                        COEF[ TabDepX[j]:TabBandNx[j]-1, *] are the ridgelet 
;                        coefficients at the scale j
;         BSIZE     --  LONG: Block size used in the ridgelet transform
;         NXB       -- LONG: Number of blocks in the x-axis direction
;         NYB       -- LONG: Number of blocks in the x-axis direction
;         OVERLAP   -- LONG: is equal to 1 if blocks are overlapping
;         HEADTRANS -- IDL STRING Array: contains the Fits Header of the decomposition
;         TABNORM   -- FLOAT Array[0:NBRSCALE-1]: Normalization value for each scale
;         TabDepX   -- Long IDL array[0:NBRSCALE-1]: Starting x-position of scale j in COEF
;         TabBandNx  -- Long IDL array[0:NBRSCALE-1]: number of columns per scale
;         TabbandNY  -- Long IDL array[0:NBRSCALE-1]: number of lines per scale
;
; KEYWORDS:
;      OPT: string which contains the differents options. Options are:
;
;         [-t type_of_ridgelet]
;              1: RectoPolar Ridgelet Transform using a standard bi-orthogonal WT 
;              2: RectoPolar Ridgelet Transform using a FFT based Pyramidal WT 
;              3: RectoPolar Ridgelet Transform using a Pyramidal WT in direct space 
;              4: Finite ridgelet transform 
;              5: Slant Stack Radon transform + FFT based pyramidal WT. 
;              6: Slant Stack Radon transformand + bi-orthogonal WT 
;              7: Slant Stack Radon transformand + pyramidal WT in direct space 
;              Default is RectoPolar Ridgelet Transform using a FFT based Pyramidal WT.
;         [-n number_of_scales]
;             Number of scales used in the ridgelet transform.
;             Default is automatically calculated. 
;         [-b BlockSize]
;             Block Size. Default is image size. 
;         [-i]
;             Print statistical information about.
;             each band. Default is no. 
;         [-O]
;             Block overlapping. Default is no. 
;         [-r]
;             Inverse RIDGELET transform.
;             NbrIter is the number of iteration used in the reconstruction of
;             the Slant Stack Radon transform. Default is 10
;         [-R NbrIter]
;             NbrIter is the number of iteration used in the reconstruction of
;             the Slant Stack Radon transform. Default is 10
;         [-x]
;             Write all bands separately as images with prefix 'band_j' (j being the band number)
;         [-z]
;             Use virtual memory.
;                default limit size: 4
;                default directory: .
;         [-Z VMSize:VMDIR]
;             Use virtual memory.
;                VMSize = limit size (megabytes) 
;                VMDIR = directory name 
;         [-v]
;             Verbose. Default is no.
;
; EXTERNAL CALLS
;           rid_trans (C++ program)
;
; EXAMPLE:
;     rid_trans, ima, Rid
;     for j=0,Rid.NBRSCALE-1 do tvscl,  Rid.coef[ Rid.TabDepX[j]:Rid.TabDepX[j]+Rid.TabBandNx[j]-1,*]
;       ridgelet transform with default options, and display all scales.
;
; HISTORY:
;       Written: Jean-Luc Starck 2005.
;       February, 2005 File creation
;-
;-----------------------------------------------------------------

function rid_getband, Rid, j
  return, Rid.coef[ Rid.TabDepX[j]:Rid.TabDepX[j]+Rid.TabBandNx[j]-1,*]
end

;-----------------------------------------------------------------

pro rid_putband, Rid, Band, j
   Rid.coef[ Rid.TabDepX[j]:Rid.TabDepX[j]+Rid.TabBandNx[j]-1,*] = Band
end

;-----------------------------------------------------------------

pro rid_trans, Ima, RidTrans, OPT=OPT 

if N_PARAMS() LT 2 then begin 
        spawn, 'rid_trans'
        print, 'CALL SEQUENCE: rid_trans, Signal, Struct_Out, OPT=Opt'
        goto, DONE
        end

Nx = (size(Ima))[1]
Ny = (size(Ima))[2]

NameIma = 'xx_signal.fits'
NameResult = 'xx_result'
NameResultFits = 'xx_result.rid'
; if TypeOrder  eq 'RING    '  or  TypeOrder  eq 'ring    ' then Order=1

writefits,  NameIma,  Ima
if not keyword_set(OPT) then OPT = " "
 
com = 'rid_trans -x  ' + OPT + ' '+ NameIma  + ' ' +  NameResult
spawn, com

Rid = readfits(NameResultFits, HeadTrans, /silent); 
;help, Rid
;help, HeadTrans
Nxrid = (size(Rid))[1]
Nyrid = (size(Rid))[2]
Nl = FXPAR(HeadTrans, "NL")
NC = FXPAR(HeadTrans, "NC")
Bsize = FXPAR( HeadTrans, "BSIZE")
Nyb = FXPAR(HeadTrans, "NLB")
Nxb = FXPAR(HeadTrans, "NCB")
Overlap = FXPAR(HeadTrans, "OVERLAP")
NbrScale = FXPAR(HeadTrans, "NBRSCALE")
TypeTrans = FXPAR(HeadTrans, "TYPE_TRA")                               

TabBandNx = lonarr(NbrScale)
TabBandNy = lonarr(NbrScale)
TabDepX = lonarr(NbrScale)

; ImaRid = fltarr(Ncrid,Nlrid)
TabDepX[NbrScale-1]=0
for j=NbrScale-1,0,-1 do begin
   NameBand = NameResult + '_band_' + STRCOMPRESS(string(j+1), /REMOVE_ALL) + '.fits'
   band = readfits(NameBand)
   delete, NameBand
   TabBandNx[j] = (size(band))[1]
   TabBandNy[j] = (size(band))[2]
   if j NE NbrScale-1 then TabDepX[j] = TabDepX[j+1] + TabBandNx[j+1]
   ; print, 'Scale ', j+1, 'Dep = ', TabDepX[j], ' Nx = ', TabBandNx[j], ' Ny = ',  TabBandNy[j]   
   ; ImaRid[TabDepX[j]:TabDepX[j]+TabBandNx[j]-1,*] = band
end

;print, ' Nl = ' , Nl, ' Nc = ', Nc, ', Bsize =', Bsize, ' NbrScale =  ' , NbrScale
;print, ' Nlb = ' , Nlb, ' Ncb = ', Ncb, ' Overlap = ', Overlap, ' TypeTrans = ', TypeTrans
;info, Rid-ImaRid
;  mr1d_struc, result, info, output, band=band

TabNorm=[0.8546,0.215, 0.1265,0.097,0.0669,0.0483,0.0342,0.0242,0.0171,0.0121]
RidTrans = {NbrScale : NbrScale, Nxrid : Nxrid, Nyrid:Nyrid, Coef : Rid, Bsize : Bsize, Nxb: Nxb, Nyb:Nyb, $
            Overlap: Overlap, HeadTrans:HeadTrans, TabDepX:TabDepX, TabBandNx: TabBandNx, TabBandNy: TabBandNy, TabNorm:TabNorm}
DONE:
 
end


pro trid, e,r, opt=opt
e = rim('einstein128.fits')
rid_trans, e, RidTrans, OPT=OPT 
hs, RidTrans

NRidTrans = RidTrans
NRidTrans.coef[*] = 0
for j=0,RidTrans.NbrScale-1 do begin
   B = rid_getband(RidTrans, j)
   rid_putband,NRidTrans, B, j
   end

info, NRidTrans.coef - RidTrans.coef 
rid_recons, RidTrans, r
end
