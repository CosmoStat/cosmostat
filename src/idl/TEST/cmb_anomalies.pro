;+
; NAME:
;        cmb_anomalies
;
; PURPOSE:
;   Computes statistical anomalies in a CMB map. 
;
; CALLING:
;     AlmANO = cmb_anomalies(Data, lmax=lmax, nsideRot=nsideRot)
;
; INPUTS:
;     Data -- IDL array of healpix map : Input image to be analyzed 
;    
; OUTPUTS:
;     AlmANO -- IDL structures with the following fields: 
;					MapName -- string : map name 
;			        lmax -- integer : lmax used in analysis
;					N -- interger: number of rotations
;					Res_LowQuadrupole -- IDL structure: results for the low quadupole analysis
;					Res_QuadOctAlignment -- IDL structure: results for the Quad-Oct alignement analysis
;					Res_PlanarOctopole -- IDL structure: results for the Oct opole Planarity analysis
;
; INPUT KEYWORDS:
;      lmax       :  integer -- lmax value in the analysis
;      nsideRot      :   resolution for angle rotations. Number of rotations is nsideRot^2*12
;
; EXTERNAL CALLS:
;       
;
; EXAMPLE:
;       Compute the anomalies for a number of rotations correcponding to nsideRot^2*12
;        AlmANO = cmb_anomalies(Data, lmax=8, nsideRot=256)
 ;         
; HISTORY:
;       Anais Rassat & Jean-Luc Starck, 2013
;       May, 2013 File creation
;--------------------------------------------------------------------------------------------------------


;======================================================================

function anomalies_lowquad, Data, nside=nside, theory=theory, norm=norm
; Written by Anais Rassat Feb 2012. 
; Adapted from code "isw_anomalies.pro Written October 2010 by Anais Rassat
; PURPOSE: Estimates quadrupole of map and returns probability that it
; is low compared to theory
;---------------------------------
; INPUT: 
; map: must be unitless (not mK or muK) overdensity
; nside: nside of input map
; Theory: theoretical value of quadrupole in l(l+1)/(2pi) muK^2
;         default is WMAP3 theory (see paper p32)
; norm: set to 1 if map is in muK
;---------------------------------
; OUTPUT: 
; prob: probability that quadrupole of input map is as low
;---------------------------------
map = Data
if keyword_set(norm) then map = map / 1d6/2.725
if not keyword_set(theory) then theory = 1252.d0 ; which is WMAP3 theory p32
if not keyword_set(nside) then nside = 32

df = 5 ; number of degrees of freedom for quadrupole = 2*l+1
delta = theory/2.d0/(2.d0+1)*2.d0*!dpi
;mrs_almtrans, map, map_alm, /tab, /complex, lmax=lmax ; alm's for ILC map

spec = mrs_powspec(map, lmax = 20)
nuk2 = 1d12*(2.725d)^2        
quad = spec[2]*nuk2 ; units = muK^2
fac = 5.d0

;prob =  imsl_chisqcdf(quad*fac/delta,df,0.d0,/double)
prob = chisqr_pdf(quad*fac/delta, df)

return, prob
end

;======================================================================

function proba_stat_t,  Nsimu,  Cl3, nsideRot
; 10000 simu and nside=128 is recommanded.
Nrot = long(nsideRot)^2 * 12L
TabStat_t = dblarr(Nrot)

; Define the Alm structure
 c = randomn(seed, nsideRot ^2*12)
mrs_almtrans, c, Alm, /complex, /tab
Alm.alm[*] = dcomplex(0,0)

t = dblarr(Nsimu)
l=3
L3 = dcomplexarr(l+1)
for i=0L, Nsimu-1 do begin
    L3[0] = dcomplex( randomn(seed) * sqrt(Cl3)  ,  0)
    for m=1,l do L3[m] = dcomplex( randomn(seed) * sqrt(Cl3/2.)  , randomn(seed) * sqrt(Cl3/2.))
    Alm.alm[l,0:l] = L3 
    mrs_almrot, Alm,  RotAlm, nsideRot=nsideRot, lmax=l,  /alm
    for r=0, Nrot-1 do begin
         L3Rot = dcomplex( RotAlm.ALMRE[l,0:l,r], RotAlm.ALMIM[l,0:l,r])
         TabL = abs(L3Rot )^2
         TabStat_t[r] =  ( 2. * TabL[3] )  /  ( 2. * total( TabL[1:*]  ) + TabL[0]  )
         t[i] = max(TabStat_t)
    end
end
      
return, t
end


;======================================================================

function axis_evil_stat,  RotAlm,  Alm2, l
    ; Alm2 = double(RotAlm.ALMre)^2 + double(RotAlm.ALMim)^2 
;   Cl0 (nˆ )	=	|a_l0 (nˆ )|^2	and C_lm(nˆ) = 2 | alm(nˆ)|^2  for m > 0 

   TabL = reform(Alm2[l,0:l, *])
   ClDen =  ( 2. * total( TabL[1:*]  ) + TabL[0]  )      
   Ratio = TabL
   Ratio[1:*, *] = 2. * Ratio[1:*, *]
   Ratio = Ratio / ClDen
   
   Rl = max( Ratio, indi )
   ind = ARRAY_INDICES(Ratio, indi)
   IndM = ind[0]
   IndPix = ind[1]
   ; PRINT, ind, Ratio[ind[0],ind[1]],  FORMAT = '(%"Value at [%d, %d] is %f")'
    ang2lb, RotAlm.THETAPHI[IndPix, 0],  RotAlm.THETAPHI[IndPix, 1] , RCoord_l,  RCoord_b     
    
   evil = {l:l, Rl: Rl, m: IndM, Coord_l: RCoord_l, Coord_b: RCoord_b, x: RotAlm.THETAPHI[IndPix, 2], y: RotAlm.THETAPHI[IndPix, 3],  z: RotAlm.THETAPHI[IndPix, 4]}
return, evil
end

;======================================================================

pro axis_evil_theta,  TabS, MeanAngle, Std
  vs = size(TabS)
  Lmax = vs[1] + 1
  Lmin = TabS[0].l

   for l= Lmin, Lmax do begin
   for l1=l+1, Lmax do begin
       S1 = TabS[l-2]
       S2 = TabS[l1-2]
       dotproduct =S1.x*S2.x + S1.y*S2.y + S1.z*S2.z
      Testangle = acos(dotproduct) * 180d / !DPI
      if Testangle GT 90 then Testangle = 180. - Testangle
      if l EQ Lmin and l1 EQ l+1 then TabTheta = Testangle $
      else TabTheta = [TabTheta, Testangle]
   endfor
   endfor
   
   
   
   MeanAngle = mean(TabTheta )
   Std = stddev(TabTheta)
end

            
 ;======================================================================

function axis_evil_proba,  Nsimu,  Spec, lmax, nsideRot
; 10000 simu and nside=128 is recommanded.
Nrot = long(nsideRot)^2 * 12L
TabStat_t = dblarr(Nrot)

; Define the Alm structure
 c = randomn(seed, nsideRot ^2*12)
mrs_almtrans, c, Alm, /complex, /tab, lmax=lmax
Alm.alm[*] = dcomplex(0,0)

TabMeanAngle = dblarr(Nsimu)
for i=0L, Nsimu-1 do begin
   ; Create the Alm of one simu
   for l=2,Lmax do begin
       Alm.alm[l,0] =  dcomplex( randomn(seed) * sqrt(Spec[l])  ,  0)
       for m=1,Lmax do  Alm.alm[l,m]  = dcomplex( randomn(seed) * sqrt(Spec[l]/2.)  , randomn(seed) * sqrt(Spec[l]/2.))
   endfor

  ; Make the rotations
   mrs_almrot, Alm,  RotAlm, nsideRot=nsideRot, lmax= lmax,  /alm
   Alm2 = double(RotAlm.ALMre)^2 + double(RotAlm.ALMim)^2 

    for l=2,lmax do begin
            S_evil =  axis_evil_stat(RotAlm,  Alm2,  l)
            if l EQ 2 then TabS = replicate( S_evil, lmax-1)
            TabS[l-2] = S_evil
     endfor
     axis_evil_theta,  TabS, MeanAngle, Std
     TabMeanAngle[i] = MeanAngle
endfor
       
return, TabMeanAngle
end


;======================================================================

pro s_parity, RotAlm,   lmax,  Spec, Smap, Sp, Sm
; reommanded nsiderot=64 and  nsideMap=8
 
   Alm2 = double(RotAlm.ALMre)^2 + double(RotAlm.ALMim)^2 
   NRot = RotAlm.NRot
   S = 0d
   Stest = dblarr(NRot)
   print, Nrot
   for l= 2, Lmax do  for m=-l,l do  Stest = Stest +   (-1d)^(l+m) * reform(Alm2[l,abs(m), *] ) / Spec[l]
   help, Stest
   Stest = Stest - (lmax-1)
   Smap = Stest
   Sp = (max(Smap) - mean(Smap)) / stddev(smap)
   Sm = abs(min(Smap) - mean(Smap)) / stddev(smap)

 end

;======================================================================

pro  s_parity_proba,  Nsimu,  Spec, lmax, nsideRot, TabSp, TabSm
Nrot = long(nsideRot)^2 * 12L
TabStat_t = dblarr(Nrot)

; Define the Alm structure
 c = randomn(seed, nsideRot ^2*12)
mrs_almtrans, c, Alm, /complex, /tab, lmax=lmax
Alm.alm[*] = dcomplex(0,0)

TabSp = dblarr(Nsimu)
TabSm = dblarr(Nsimu)
for i=0L, Nsimu-1 do begin
   ; Create the Alm of one simu
   for l=2,Lmax do begin
       Alm.alm[l,0] =  dcomplex( randomn(seed) * sqrt(Spec[l])  ,  0)
       for m=1,Lmax do  Alm.alm[l,m]  = dcomplex( randomn(seed) * sqrt(Spec[l]/2.)  , randomn(seed) * sqrt(Spec[l]/2.))
   endfor

  ; Make the rotations
   mrs_almrot, Alm,  RotAlm, nsideRot=nsideRot, lmax= lmax,  /alm
   s_parity, RotAlm,   lmax,  Spec, Smap, Sp, Sm
   TabSp[i] = Sp
   TabSm[i] = Sm
endfor
 end
          
;======================================================================
; Main routines
;======================================================================

function cmb_anomalies, Data, lmax=lmax, nsideRot=nsideRot, TheoryQuadrupole=TheoryQuadrupole, norm=norm, verb=verb, AnoQuadThe=AnoQuadThe,  AnoQuadOctAlignment=AnoQuadOctAlignment,  AnoPlanarOctopole=AnoPlanarOctopole, RotAlm=RotAlm, MapName=MapName
; Written by Anais Rassat & J.L Starck Oct 2012. 
;  
;---------------------------------
; INPUT: 
; map: must be unitless (not mK or muK) overdensity
; Theory: theoretical value of quadrupole in l(l+1)/(2pi) muK^2
;         default is WMAP3 theory (see paper p32)
; norm: set to 1 if map is in muK
;---------------------------------
; OUTPUT: 
; prob: probability that quadrupole of input map is as low
;---------------------------------

map = Data
N = N_elements(map)

if not keyword_set(lmax) then lmax = 5L
if keyword_set(Verb) then print, 'Lmax = ',  lmax
; if keyword_set(norm) then map = map / 1d6/2.725
if not keyword_set(TheoryQuadrupole) then TheoryQuadrupole = 1252.d0 ; which is WMAP3 theory p32
; if not keyword_set(nsideRot) then nsideRot = 128L
if not keyword_set(nsideRot) then nsideRot = 32L
N =  nsideRot^2*12L

; if nsideRot ne sqrt(N/12.) then begin
;   print, 'cmb_anomalies: check size of map with nsideRot'
;stop
; endif

 Spec = mrs_powspec(map, lmax=20)
if not keyword_set(RotAlm) then  mrs_almrot, map, RotAlm, lmax=lmax, nsideRot=nsideRot
 Alm2 = double(RotAlm.ALMre)^2 + double(RotAlm.ALMim)^2 
; Alm2 = RotAlm.ALMre^2 + RotAlm.ALMim^2 
if keyword_set(verb) then print, 'Calculated rotations up to lmax = ', lmax

;=======================   The low quadrupole ======================
print, '-------------------------------------'
print, 'BEGINNING: calculation of low quadrupole anomaly'
print, '-------------------------------------'
AnoQuadThe=1

theory = TheoryQuadrupole
DataQuadrupole = 0.
ProbLowQuadrupole = 0.
if keyword_set(AnoQuadThe) then begin 
   df = 5 ; number of degrees of freedom for quadrupole = 2*l+1
   delta = theory/2.d0/(2.d0+1)*2.d0*!dpi
   nuk2 = 1d12*(2.725d)^2        
   DataQuadrupole = spec[2]*nuk2 ; units = muK^2
   fac = 5.d0
   ;prob =  imsl_chisqcdf(quad*fac/delta,df,0.d0,/double)
   ProbLowQuadrupole = chisqr_pdf(DataQuadrupole*fac/delta, df)
   if keyword_set(Verb) then print, 'Prob(DataQuad | TheoQuad) = ', ProbLowQuadrupole
   
   Res_LowQuadrupole= {TheoryQuadrupole: TheoryQuadrupole,   DataQuadrupole: DataQuadrupole,  ProbLowQuadrupole: ProbLowQuadrupole}
end


;=======================  The quadrupole-octopole alignment ======================
print, '-------------------------------------'
print, 'BEGINNING: calculation of quadrupole-octopole anomaly'
print, '-------------------------------------'


AnoQuadOctAlignment = 1
Verb=1
TabQl=1

QlMax = {Ind: 0L,  Val: 0d, Theta: 0d, Phi: 0d, l: 0d, b:0d, x:0d, y:0d, z:0d}

; Find the preferred axis for l=2 and l=3
;  to speed we can put nsideRot < 512, nsideRot=512 gives a precision of < 1' for the preferred axis (check in Oliveira paper)
if  keyword_set(AnoQuadOctAlignment) then begin 
    TabSum = dblarr(N)
    TabQlMax = replicate(QlMax, lmax+1)
    for l=2,lmax do begin ; calculate statistic for each l
      TabL = reform(Alm2[l,0:l, *]) ; |alm|^2 for a given l and for pixel i (last variable
      indm = lindgen(l+1)
      indm2 = double(indm * indm) ; m^2
      for i=0L,N-1 do TabSum[i] = 2.d0 * total( TabL[*,i] * indm2 ); for each pixel value, calculate stat, added double precision
      
      TabQlMax[l].Val = max(TabSum, indAngle)
      TabQlMax[l].Ind = indAngle

      TabQlMax[l].Theta  = RotAlm.THETAPHI[indAngle, 0] ; angle theta in radian
      TabQlMax[l].Phi  = RotAlm.THETAPHI[indAngle, 1] ; angle phi in radia
      ang2lb, TabQlMax[l].Theta ,  TabQlMax[l].Phi, Coord_l,  Coord_b
      TabQlMax[l].l  = Coord_l
      TabQlMax[l].b  = Coord_b
      TabQlMax[l].x = RotAlm.THETAPHI[indAngle, 2]
      TabQlMax[l].y  = RotAlm.THETAPHI[indAngle, 3]
     TabQlMax[l].z  = RotAlm.THETAPHI[indAngle, 4]
   endfor
     DotProduct =  TabQlMax[2].x *  TabQlMax[3].x + TabQlMax[2].y *  TabQlMax[3].y  + TabQlMax[2].z *  TabQlMax[3].z
     Angle = acos(DotProduct)*180d / !DPI
     Proba = (1d - DotProduct)*100d
     if keyword_set(Verb) then print,  'n_2 . n_3 = ', DotProduct,  ', Angle (deg) = ', Angle, ', Probability (%) = ', Proba
   ; write(*,*)  dotproduct, acos(dotproduct)*rad2deg, (1.d0 - dotproduct)*100.d0
    ; if keyword_set(Verb) then for l=2,lmax do  print, 'l= ', strc(l), ' Ql = ', strc(TabQl[l,0]), ', Angle (l,b) = (', strc(TabQl[l,3]) , ', ', strc(TabQl[l,4]) , ')'
    Res_QuadOctAlignment = { TabQlMax: TabQlMax, DotProduct: DotProduct , Angle : Angle , Proba: Proba}
endif

;=======================  The planar octopole ======================
print, '-------------------------------------'
print, 'BEGINNING: calculation of octopole planarity anomaly'
print, '-------------------------------------'
AnoPlanarOctopole=1
Nsimu=1000L
 Cl3 = Spec[3]
 
if  keyword_set(AnoPlanarOctopole) then begin 
    l=3
    TabL = reform(Alm2[l,0:l, *])
    TabM3 = dblarr(N)
    for i=1L,N-1 do TabM3[i] =  ( 2. * TabL[3,i] )  /  ( 2. * total( TabL[1:*,i]  ) + TabL[0,i]  )
    StatM3 = max(TabM3, indAngle)
    IndM3 = indAngle
    M3Theta  = RotAlm.THETAPHI[IndM3, 0]
    M3Phi  = RotAlm.THETAPHI[IndM3, 1]
    M3x = RotAlm.THETAPHI[indAngle, 2]
    M3y  = RotAlm.THETAPHI[indAngle, 3]
    M3z  = RotAlm.THETAPHI[indAngle, 4]
    ang2lb, M3Theta,  M3Phi , M3Coord_l,  M3Coord_b
    
    TabSimT = proba_stat_t(Nsimu,  Cl3, nsideRot)
    ; 10000 simu and nside=128 is recommanded.

    ind = where(  StatM3  LT TabSimT,  count)
   if count EQ 0 then print, "PB in calculation of octopole planarity anomaly … "
   Proba = (1. - double(count) /double(Nsimu)) * 100.
   Res_PlanarOctopole = { Stat_t: StatM3,   IndAngle: IndM3, Theta: M3Theta,  Phi : M3Phi ,  Coord_l : M3Coord_l , Coord_b : M3Coord_b , x: M3x,  y: M3y, z: M3z, Proba: Proba}
   if keyword_set(Verb) then print,  'Stat_t = ', Res_PlanarOctopole.Stat_t,  ',  Probability (%) = ', Res_PlanarOctopole.Proba
 end

;=======================  Evil Axis  ======================
print, '-------------------------------------'
print, 'BEGINNING: calculation of ''axis of evil'' anomaly'
print, '-------------------------------------'
AnoAxisOfEvil=1
Res_AvilAxis=0
 if  keyword_set(AnoAxisOfEvil) then begin 
     for l=2,lmax do begin
            S_evil =  axis_evil_stat(RotAlm,  Alm2,  l)
            if l EQ 2 then TabS = replicate( S_evil, lmax-1)
            TabS[l-2] = S_evil
     endfor

    axis_evil_theta,  TabS, MeanAngle, Std
 
   Proba = 0.
   ; Nsimu=1000L
   Nsimu = 5L;
   TabMeanAngle = axis_evil_proba(Nsimu,  Spec, lmax, nsideRot)
   ind = where(  MeanAngle  LT TabMeanAngle,  count)
   if count EQ 0 then print, "PB in calculation of axis of evil anomaly … "
   Proba = (double(count) /double(Nsimu)) * 100.
   if keyword_set(Verb) then print,  ' Angle: ', MeanAngle, ' +/- ', Std, ', Proba = ', Proba
   Res_AOE = {  TabS: TabS,   $
                                MeanAngle: MeanAngle, $
                                Std: Std, $
                                TabMeanAngle: TabMeanAngle, $
                                Proba: Proba}
endif

;=======================  Mirror Parity  ======================

AnoSParity=1
Res_SParity = 0

NsideInMapSParity = 8L
NsideRotSParity = 64L
Nsimu = 5L;

MapLowResol = mrs_resize(map, nside=NsideInMapSParity)
Spec = mrs_powspec(MapLowResol, lmax=20)
mrs_almrot, MapLowResol, RotAlm, lmax=lmax, nsideRot=NsideRotSParity
 
 if  keyword_set(AnoSParity) then begin 
      s_parity, RotAlm,  lmax,  Spec,  Smap, Sp, Sm
      Proba = 0.
      s_parity_proba,  Nsimu,  Spec, lmax, nsideRot, TabSp, TabSm
      ind = where(  Sp  LT TabSp,  count)
      if  count EQ 0 then print, "PB in calculation of Mirror Parity anomaly … "
      ProbaSp = (1. - double(count) /double(Nsimu)) * 100.

      ind = where(  Sm  LT TabSm,  count)
      if  count EQ 0 then print, "PB in calculation of Mirror Parity anomaly … "
      ProbaSm = (1. - double(count) /double(Nsimu)) * 100.
      
    if keyword_set(Verb) then print,  ' Sp = ', Sp,  ', Proba = ', ProbaSp
    if keyword_set(Verb) then print,  ' Sm = ', Sm,  ', Proba = ', ProbaSm


      Res_MirrorParity = {Smap: Smap, $
                                  Sp: Sp, $
                                  Sm: Sm, $
                                  TabSp: TabSp, $
                                  TabSm: TabSm, $
                                ProbaSp: ProbaSp, $
                                ProbaSm: ProbaSm}       
endif


;=======================  The S-Stat  ======================
if not keyword_set(MapName) then MapName = 'Map'
Anomalies = { MapName: MapName, $
	                 lmax:lmax, N:N,  $
                     Res_LowQuadrupole: Res_LowQuadrupole, $  ; Quadrupole
                     Res_QuadOctAlignment: Res_QuadOctAlignment, $  ; Quad-Oct alignement
                     Res_PlanarOctopole: Res_PlanarOctopole, $         ; Oct Planarity
                     Res_AOE: Res_AOE, $         ; Axis of evil
                     Res_MirrorParity: Res_MirrorParity $  ; Mirror Parity 
                   }

return, Anomalies
end
