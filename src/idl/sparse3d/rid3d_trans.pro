;+
; NAME: 
;       RID3D_TRANS
;
; PURPOSE:
;       Compute the 3D ridgelet  transform of a cube. The output is an IDL structure.
;       A band at scale j (j=0..NBRSCALE-1) can be extracted using the function
;       function rid3d_getband(Rid, j) (ex: Scale2 = rid3d_getband(RidTrans, 2))
;       and a band can be inserted in the transformation using the routine  rid3d_putband
;       (ex:  rid3d_putband, RidTrans, Scale2, 2).
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
;         NXRID     -- LONG: Number of pixels in the x-axis direction
;         NYRID     -- LONG: Number of pixels in the y-axis direction
;         COEF      -- 2D IDL array: Ridgelet coefficients
;                             Image which contains all ridgelet coefficients
;                        COEF[ TabDepX[j]:TabBandNx[j]-1, *] are the ridgelet 
;                        coefficients at the scale j
;         BSIZE     --  LONG: Block size used in the ridgelet transform
;         NXB       -- LONG: Number of blocks in the x-axis direction
;         NYB       -- LONG: Number of blocks in the y-axis direction
;         NZB       -- LONG: Number of blocks in the z-axis direction
;         OVERLAP   -- LONG: is equal to 1 if blocks are overlapping
;         HEADTRANS -- IDL STRING Array: contains the Fits Header of the decomposition
;         TABNORM   -- FLOAT Array[0:NBRSCALE-1]: Normalization value for each scale
;         TabDepX   -- Long IDL array[0:NBRSCALE-1]: Starting x-position of scale j in COEF
;         TabBandNx  -- Long IDL array[0:NBRSCALE-1]: number of columns per scale
;         TabbandNY  -- Long IDL array[0:NBRSCALE-1]: number of lines per scale
;
; KEYWORDS:
;      OPT: string which contains the differents options. Options are:
;         [-t type_of_ridgelet]
;              1: Standard bi-orthogonal WT 
;              2: FFT based Pyramidal WT 
;              Default is FFT based Pyramidal WT.
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
;         [-v]
;             Verbose. Default is no.
;
; EXTERNAL CALLS
;           rid3d_trans (C++ program)
;
; EXAMPLE:
;     rid3d_trans, Cube, Rid
;     for j=0,Rid.NBRSCALE-1 do tvscl,  Rid.coef[ Rid.TabDepX[j]:Rid.TabDepX[j]+Rid.TabBandNx[j]-1,*]
;       ridgelet transform with default options, and display all scales.
;
; HISTORY:
;       Written: Jean-Luc Starck 2005.
;       September, 2005 File creation
;-
;-----------------------------------------------------------------

function rid3d_getband, Rid, j
  return, Rid.coef[ Rid.TabDepX[j]:Rid.TabDepX[j]+Rid.TabBandNx[j]-1,*]
end

;-----------------------------------------------------------------

pro rid3d_putband, Rid, Band, j
   Rid.coef[ Rid.TabDepX[j]:Rid.TabDepX[j]+Rid.TabBandNx[j]-1,*] = Band
end

;-----------------------------------------------------------------

pro rid3d_trans, Ima, RidTrans, OPT=OPT 

if N_PARAMS() LT 2 then begin 
        spawn, 'rid3d_trans'
        print, 'CALL SEQUENCE: rid3d_trans, Signal, Struct_Out, OPT=Opt'
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
 
com = 'rid3d_trans -x  ' + OPT + ' '+ NameIma  + ' ' +  NameResult
spawn, com

Rid = readfits(NameResultFits, HeadTrans, /silent); 
;help, Rid
;help, HeadTrans
Nxrid = (size(Rid))[1]
Nyrid = (size(Rid))[2]
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
RidTrans = {NbrScale : NbrScale, Nxrid : Nxrid, Nyrid:Nyrid, Coef : Rid, Bsize : Bsize, Nxb: Nxb, Nyb:Nyb, Nzb: Nzb, $
            Overlap: Overlap, HeadTrans:HeadTrans, TabDepX:TabDepX, TabBandNx: TabBandNx, TabBandNy: TabBandNy, TabNorm:TabNorm}
DONE:
 
end

pro run_simu, RidMean, RidSigma, bs=bs, dens=dens, n=n, pois=pois,gauss=gauss
; run_simu, RidMean, RidSigma, dens=3., n=10.,/pois, bs=8 
; run_simu, RidMean, RidSigma,  /pois, bs=8 

if not keyword_set(bs) then bs=8
if not keyword_set(n) then n=1000
if not keyword_set(dens) then dens= double(3870.)/ 1093440.
print, 'Block size = ', bs, ' Dens = ', dens
seed=long([10])
for i=0,n-1 do begin
   if i mod 20 eq 0 then print, 'SIMU ', i
   cube = fltarr(bs,bs,bs)
   cube[*] = dens
   ; help, seed, seed[0]
   if keyword_set(pois) then cube = poisson_image(cube, seed[0]) $
   else if keyword_set(gauss) then cube = randomn(seed, bs,bs,bs) * dens
   rid3d_trans, cube, RidTrans, OPT=OPT 
   Coef = double(RidTrans.coef)
   if i EQ 0 then begin
        RidMean =  Coef
	RidSigma = Coef*Coef
   end else begin
       RidMean = RidMean + Coef
       RidSigma = RidSigma + Coef*Coef
   end
end
RidMean = RidMean / float(n)
RidSigma = RidSigma / float(n)
RidSigma =  RidSigma - RidMean*RidMean
ind = where(RidSigma LT 0, c)
if c GT 0 then RidSigma(ind) = 0
RidSigma =  sqrt(RidSigma)

vs = size(RidSigma)
RidNorm = fltarr(vs[1],vs[2],2)
RidNorm[*,*,0] = RidMean
RidNorm[*,*,1] = RidSigma
Name='rid3dnorm'+STRCOMPRESS(string(bs), /REMOVE_ALL) + '.fits'
writefits, Name, RidNorm
end

pro trid3d, e,r, opt=opt
N=8
e = fltarr(N,N,N)
e(*,N/2,*)=1
e = readfits('b.fits')
rid3d_trans, e, RidTrans, OPT=OPT 
 
NRidTrans = RidTrans
NRidTrans.coef[*] = 0
for j=0,RidTrans.NbrScale-1 do begin
   B = rid3d_getband(RidTrans, j)
   rid3d_putband,NRidTrans, B, j
   end

info, NRidTrans.coef - RidTrans.coef 
rid3d_rec, RidTrans, r
end
