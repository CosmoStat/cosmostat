function mrs_hc_wtima, Data, NbrScale=NbrScale, HCCoef=HCCoef, HCStat=HCStat, NoiseMap=NoiseMap, Cl=Cl, HCMax=HCMax, MadHC=MadHC, CovMat=CovMat, WT=WT, Mask=Mask

if not keyword_set(WT) then mrs_wttrans, Data, WT,  NbrScale=NbrScale
if keyword_set(NoiseMap ) then begin
  if not keyword_set(CovMat) then CovMat = mrs_GetWTLocalVariance(NoiseMap,NbrScale=NbrScale, BS=BS,Cl=Cl) 
end

HCCoef = WT.coef
HCCoef[*]=0.
HCStat = fltarr(WT.NBRSCALE)
MadHC = fltarr(WT.NBRSCALE) 

for j=0, WT.NBRSCALE-2 do begin
   Scale = WT.coef[*,j]
   Np = float(N_ELEMENTS(Scale))

   if keyword_set(Mask) then begin
      IndMask = where(Mask EQ 1,  Np)
      Scale = Scale[indMask]
      Np = float(Np)
   end
   if keyword_set(NoiseMap) then begin
      VarData =  total(Scale^2)
      CovMatScale = CovMat[*,j]
      if keyword_set(Mask) then CovMatScale = CovMatScale[indMask]
      VarNoise = total(CovMatScale)
      VarCMB = (VarData - VarNoise) / Np
      if VarCMB LT 0 then VarCMB = 0
      Scale = WT.coef[*,j] / sqrt(CovMat[*,j] + VarCMB)
      
      HCT, Scale, HCI, HCI2, /zeromean, HCIma=HCIma, /nonorm
      ; TabStat =  get_stat(Scale, HCIma=HCIma)
      ; HCI2  = TabStat[6]
   end else  begin
      TabStat =  get_stat(Scale, HCIma=HCIma)
      HCI2  = TabStat[6]
   end
   MadHC[j] = mad( HCIma )
   HCIma = alog(1.+abs(HCIma))
   HCStat[j] = HCI2
   HCCoef[*,j] = HCIma   ;   /  MadHC[j]
   tvs, HCCoef[*,j], tit='Scale ' + strc(j+1), win=1
   print, 'Scale ', j+1, ', Mad = ',  MadHC[j]

   if j EQ 0 then HCMax = HCCoef[*,j] $
   else begin
     ind = where ( abs(HCIma) GT abs(HCMax), c)
     if c GT 0 then HCMax[ind] = HCIma[ind]
   end
end
W = WT
W.coef = abs(HCCoef)
mrs_wtrec, w, HCRec

return, HCRec / float(W.NBRSCALE-1)
end

;=====================================================================================

function planck_gen_noise, WT_CovMat=WT_CovMat, NoiseMap=NoiseMap, NbrScale=NbrScale, Cl=Cl, BS=BS, Niter=Niter
  
if not keyword_set(Niter) then Niter = 10
if not keyword_set(WT_CovMat) then begin
if not keyword_set(NoiseMap) then begin
    print, 'CALLING SEQUENCE: NoiseRea = planck_gen_noise(WT_CovMat=WT_CovMat, NoiseMap=NoiseMap, NbrScale=NbrScale)'
    goto, DONE
    end  else begin
    print, "OK"
     WT_CovMat = mrs_GetWTLocalVariance(NoiseMap,NbrScale=NbrScale, BS=BS,Cl=Cl, WTTrans=WTTrans)
   end
end

NoiseRea = randomn(seed, WTTrans.npix)
 Niter =10
for i=0, Niter-1 do begin
      WT_ReaCovMat = mrs_GetWTLocalVariance(NoiseRea,NbrScale=NbrScale, BS=BS,Cl=Cl, WTTrans=WTTrans)
      for j=0, WTTrans.NBRSCALE-1 do begin
           Scale = WTTrans.coef[*,j]
           WTTrans.coef[*,j] = Scale /  sqrt(WT_ReaCovMat[*,j])  * sqrt(WT_CovMat[*,j]) 
    end
    mrs_wtrec, WTTrans, NoiseRea
end


DONE:

return, NoiseRea

end

;=====================================================================================

