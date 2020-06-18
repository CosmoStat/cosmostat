;+
; NAME: 
;       MRS_RIDREC
;
; PURPOSE: 
;        Reconstruct an image on the Sphere from its ridgelet transform (see mrs_ridtrans).   
;
; CALLING:
;       MRS_RIDREC, Rid_Struct, result
;
; INPUT:
;       Rid_Struct : IDL structure; Ridgelet transform structure (see MRS_RIDTRANS) 
;          
; OUTPUTS:
;      Result:  1D array of an Healpix image (NESTED format)
;             
; EXTERNAL CALLS
;           invrid2d
;
; EXAMPLE:
;   mrs_ridtrans, Imag, Rid 
;   mrs_ridrec, Rid, RecIma
;       Ridgelet transform and reconstruction
;
; HISTORY:
;       Written:  Jean-Luc Starck & Yassir Moudden & Ludovic Poupard, 2005.
;       February, 2005 File creation
;-
;-----------------------------------------------------------------

pro mrs_ridrec, Trans, HealpixIma 

if N_PARAMS() LT 2 then begin 
        spawn, 'rid_trans'
        print, 'CALL SEQUENCE: mrs_ridrec, Rid_Struct, result'
        goto, DONE
        end
        
NameSig = 'xx_signal.fits'
NameRid = 'xx_coef.rid'
 
for f=0,11 do begin
  if Trans.bin EQ 1 then BEGIN
    Coef = Trans.Coef[*,*,f]
    RidTrans = {NbrScale : Trans.NbrScale, Nxrid : Trans.Nxrid, Nyrid:Trans.Nyrid, Coef : Coef, Bsize : Trans.Bsize, Nxb: Trans.Nxb, Nyb:Trans.Nyb, $
            Overlap: Trans.Overlap, HeadTrans:Trans.HeadTrans, TabDepX:Trans.TabDepX, TabBandNx: Trans.TabBandNx, TabBandNy: Trans.TabBandNy, TabNorm:Trans.TabNorm}
    writefits,  NameRid, Coef, Trans.HEADTRANS
    com = 'rid_trans -r ' + NameRid  + ' ' +  NameSig
    spawn, com
    Rec = readfits(NameSig, /silent); 
    END ELSE BEGIN
       RidTrans = Trans.Coef[f]
       invrid2d, RidTrans, Rec 
    END
    if f eq 0 then   CubeFace = fltarr((size(Rec))[1], (size(Rec))[1], 12)
    CubeFace[*,*,f] = Rec 
end

put_all_faces, CubeFace, HealpixIma

delete, NameSig
delete, NameRid 
DONE:
end


pro fig_backproj, RidTrans, b
; Example of backproject of Dirac
i = fltarr(196608L)
mrs_ridtrans, i, RidTrans, opt='-b 128 -O -t5'
RidTrans.coef[*] = 0

NS = RidTrans.NBRSCALE-1 
print, ' NBRSCALE = ', NS
for j=0,NS-1 do begin
  CubeBand = mrs_ridget(RidTrans, j) 
  vs = size(CubeBand)
  Nx = vs[1] / RidTrans.nxb
  Ny = vs[2] / RidTrans.nyb
  CubeBand[Nx/2, (Ny-10)/(NS-j), j] = 100. / RidTrans.TabNorm[j]
  
  CubeBand[Nx/2, (Ny-20)/(j+1), NS+j] = 100. / RidTrans.TabNorm[j]
  
  mrs_ridput, RidTrans, CubeBand, j
  end
 
mrs_ridrec,  RidTrans, b
tvs, b
end
