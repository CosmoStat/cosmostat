;+
; NAME:
;        stat_anomalies
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

function stat_parity,map,NsideInMapSParity=NsideInMapSParity, NsideRotSParity=NsideRotSParity,nsimu = nsimu,lmax=lmax, noprob=noprob, noSmap=noSmap
; Written by Anais Rassat & Jean-Luc Starck, May 2013
; Adapted from code anomalies_parity.f90 written July 2012 by Anais Rassat
; PURPOSE: Estimates odd and even mirror parity parameters S+ and S-,
; as well as parity map S(\hat{n}) as well as probability of this occuring by chance
;---------------------------------
; INPUT: 
;---------------------------------
; OUTPUT: 
;---------------------------------
; CHECKED: Reproduces similar results to Rassat & Starck for W3 ILC
; map, with about 1% accuracy. Not sure yet where this difference
; comes from.
;---------------------------------
  if not keyword_set(lmax) then lmax = 5L
  if not keyword_set(nsimu) then nsimu = 100L
  if not keyword_set(NsideInMapSParity) then NsideInMapSParity = 8L
  if not keyword_set(NsideRotSParity) then NsideRotSParity = 64L
  
  MapLowResol = mrs_resize(map, nside=NsideInMapSParity)
  Spec = mrs_powspec(MapLowResol, lmax=20)
  mrs_almrot, MapLowResol, RotAlm, lmax=lmax, nsideRot=NsideRotSParity
  
  s_parity, RotAlm,  lmax,  Spec,  Smap, Sp, Sm
  if not keyword_set(noprob) then begin
     Proba = 0.
     s_parity_proba,  Nsimu,  Spec, lmax, NsideRotSParity, TabSp, TabSm
     ind = where(  Sp  LT TabSp,  count)
     if  count EQ 0 then print, "stat_parity: PB in calculation of Mirror Parity anomaly … "
     ProbaSp = double(count)/double(Nsimu)*100. ; corrected by Anais Rassat 21.3.2014
;     ProbaSp = (1.d - double(count) /double(Nsimu)) * 100.     ; corrected by Anais Rassat 20/6/13
     ind = where(  Sm  LT TabSm,  count)
     if  count EQ 0 then print, "stat_parity: PB in calculation of Mirror Parity anomaly … "
     ProbaSm = double(count)/double(Nsimu)*100. ; corrected by Anais Rassat 21.3.2014
;     ProbaSm = (1.d - double(count) /double(Nsimu)) * 100.     ; corrected by Anais Rassat 20/6/13 
 endif else begin
     ProbaSp = -1.
     ProbaSm = -1.
     TabSp = -1.
     TabSm = -1.
  endelse
  
  if keyword_set(noSmap) then Smap=0
  
  Res_MirrorParity = {Smap: Smap, $
                      Sp: Sp, $
                      Sm: Sm, $
                      TabSp: TabSp, $
                      TabSm: TabSm, $
                      ProbaSp: ProbaSp, $
                      ProbaSm: ProbaSm}       
  return, Res_MirrorParity
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
    Rcoord_l = 180.d - Rcoord_l ; to use same units as de Oliveira-Costa & Land/Maguiejo
    
   evil = {l:l, Rl: Rl, m: IndM, Coord_l: RCoord_l, Coord_b: RCoord_b, x: RotAlm.THETAPHI[IndPix, 2], y: RotAlm.THETAPHI[IndPix, 3],  z: RotAlm.THETAPHI[IndPix, 4]}
return, evil
end
;======================================================================
function axis_evil_proba,  Spec, lmax, nsideRot,nsimu=nsimu
if not keyword_set(nsimu) then nsimu = 100L
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
function stat_aoe,Spec,Alm2,RotAlm,lmax,nsimu=nsimu,noprob=noprob
; Written by Anais Rassat & Jean-Luc Starck, May 2013
; Adapted from code anomalies_aoe.f90 written June 2012 by Anais Rassat
; PURPOSE: Estimates `Axis of Evil' for l=2, lmax and probability of this occuring by chance
;---------------------------------
; INPUT: 
;---------------------------------
; OUTPUT: 
;---------------------------------
; CHECKED: Reproduces similar results both in preferred phase and
; direction as in Land & Magueijo 2003 for W3 ILC map (without
; subtraction of kinetic Doppler quadrupole)
;---------------------------------
if not keyword_set(nsimu) then   Nsimu = 100L
nsiderot = RotAlm.nsiderot
 for l=2,lmax do begin
      S_evil =  axis_evil_stat(RotAlm,  Alm2,  l)
      if l EQ 2 then TabS = replicate( S_evil, lmax-1)
      TabS[l-2] = S_evil
   endfor
   
   axis_evil_theta,  TabS, MeanAngle, Std
   
   if not keyword_set(noprob) then begin
      Proba = 0.
      TabMeanAngle = axis_evil_proba(Spec, lmax, nsideRot,Nsimu=Nsimu)
      
      ind = where(  MeanAngle  LT TabMeanAngle,  count)
      if count EQ 0 then print, "stat_aoe: PB in calculation of axis of evil anomaly … "
      Proba = (1. - double(count) /double(Nsimu)) * 100. ;corrected by Anais Rassat 20/6/13
   endif else begin
      proba = -1.
      tabmeanangle = -1.
   endelse
   Res_AOE = {  TabS: TabS,   $
                MeanAngle: MeanAngle, $
                Std: Std, $
                TabMeanAngle: TabMeanAngle, $
                Proba: Proba}
return, Res_AOE
end

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
     print, i
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
function stat_planar_oct, spec, alm2, rotalm, nsimu=nsimu,noprob=noprob
; Written by Anais Rassat & Jean-Luc Starck, May 2013
; Adapted from code anomalies_octplan.f90 written June 2011 by Anais Rassat
; PURPOSE: Estimates planarity of l=3 and probability of this occuring by chance
;---------------------------------
; INPUT: 
;---------------------------------
; OUTPUT: 
;---------------------------------
; CHECKED: Reproduces same results as Rassat, Starck & Dupé 2013 for
; W3 ILC map
;---------------------------------
  if not keyword_set(nsimu) then nsimu = 4L
  nsiderot = long(sqrt(RotAlm.nrot/12L))
  Cl3 = Spec[3]
  l=3
  TabL = reform(Alm2[l,0:l, *])
  TabM3 = dblarr(RotAlm.nrot)
  for i=1L,RotAlm.nrot-1 do TabM3[i] =  ( 2. * TabL[3,i] )  /  ( 2. * total( TabL[1:*,i]  ) + TabL[0,i]  )
  StatM3 = max(TabM3, indAngle)
  IndM3 = indAngle
  M3Theta  = RotAlm.THETAPHI[IndM3, 0]
  M3Phi  = RotAlm.THETAPHI[IndM3, 1]
  M3x = RotAlm.THETAPHI[indAngle, 2]
  M3y  = RotAlm.THETAPHI[indAngle, 3]
  M3z  = RotAlm.THETAPHI[indAngle, 4]
  ang2lb, M3Theta,  M3Phi , M3Coord_l,  M3Coord_b
  if not keyword_set(noprob) then begin
;     print, 'stat_planar_oct: calculating probability using
;     simulations.'
     stop
     TabSimT = proba_stat_t(Nsimu,  Cl3, nsideRot)
                                ; 10000 simu and nside=128 is recommanded.
     ind = where(  StatM3  GT TabSimT,  count)
     if count EQ 0 then print, "stat_planar_oct: PB in calculation of octopole planarity anomaly … "
     Proba = ( 1. - double(count) /double(Nsimu)) * 100. 
  endif else begin
     proba = -1.
     tabsimT = -1.
;     print, 'stat_planar_oct: not calculating probability of t statistic.'
  endelse
  Res_PlanarOctopole = { Stat_t: StatM3,   IndAngle: IndM3, Theta: M3Theta,  Phi : M3Phi ,  Coord_l : M3Coord_l , Coord_b : M3Coord_b , x: M3x,  y: M3y, z: M3z, Proba: Proba, t_simu:TabsimT}
return, Res_PlanarOctopole
end

;======================================================================
function stat_l2l3, map, alm2, RotAlm
; Written by Anais Rassat & Jean-Luc Starck, May 2013
; Adapted from code anomalies_l2l3_axis.f90 written Feb 2011 by Anais Rassat
; PURPOSE: Estimates preferred direction for l=2 and l=3 and the angle
; between the two axes, also returns probability of this occuring by chance
;---------------------------------
; INPUT: 
;---------------------------------
; OUTPUT: 
;---------------------------------
; CHECKED: Reproduces same results as Rassat, Starck & Dupé 2013 for
; W3 ILC map
;---------------------------------
  lmax = 3
  TabQl=1
  QlMax = {l:0L,Ind: 0L,  Val: 0d, Theta: 0d, Phi: 0d, Coord_l: 0d, Coord_b:0d, x:0d, y:0d, z:0d}
  TabSum = dblarr(RotAlm.Nrot)
  TabQlMax = replicate(QlMax, lmax+1)
  for l=2,lmax do begin             ; calculate statistic for each l
     TabL = 0.d0
     TabL = reform(Alm2[l,0:l, *])  ; |alm|^2 for a given l and for pixel i (last variable
     indm = lindgen(l+1)
     indm2 = double(indm * indm)                                    ; m^2
     for i=0L,RotAlm.Nrot-1 do TabSum[i] = 2.d0 * total( TabL[*,i] *indm2 ) ; for each pixel value, calculate stat, added double precision
     TabQlMax[l].l = l
     TabQlMax[l].Val = max(TabSum, indAngle)
     TabQlMax[l].Ind = indAngle
     
     TabQlMax[l].Theta  = RotAlm.THETAPHI[indAngle, 0] ; angle theta in radian
     TabQlMax[l].Phi  = RotAlm.THETAPHI[indAngle, 1]   ; angle phi in radia
     ang2lb, TabQlMax[l].Theta ,  TabQlMax[l].Phi, Coord_l,  Coord_b
     TabQlMax[l].Coord_l  = 180.d -Coord_l
     TabQlMax[l].Coord_b  = Coord_b
     TabQlMax[l].x = RotAlm.THETAPHI[indAngle, 2]
     TabQlMax[l].y  = RotAlm.THETAPHI[indAngle, 3]
     TabQlMax[l].z  = RotAlm.THETAPHI[indAngle, 4]
  endfor
  DotProduct =  abs(TabQlMax[2].x *  TabQlMax[3].x + TabQlMax[2].y *  TabQlMax[3].y  + TabQlMax[2].z *  TabQlMax[3].z)
  Angle = acos(DotProduct)*180d / !DPI
  if Angle gt 90.d0 then Angle = 180.d - Angle
  Proba = (1d - DotProduct)*100d
  Res_QuadOctAlignment = { TabQlMax: TabQlMax, DotProduct: DotProduct , Angle : Angle , Proba: Proba, SimuProb:0.}
  return, Res_QuadOctAlignment
end

;================================================================================
function stat_lowquad, spec,theoryquad=theoryquad,norm=norm,dataquad=dataquad
; Written by Anais Rassat Feb 2012. 
; Adapted from code "isw_anomalies.pro Written October 2010 by Anais Rassat
; PURPOSE: Estimates quadrupole of map and returns probability that it
; is low compared to theory
;---------------------------------
; INPUT: 
; spec: must be in muK^2
; nside: nside of input map
; Theory: theoretical value of quadrupole in muK^2
;         default is WMAP9 theory
; norm: set to 1 if map is in mK -> 
;---------------------------------
; OUTPUT: 
; prob: probability that quadrupole of input map is as low
;---------------------------------
; CHECKED: Reproduces same results as Rassat, Starck & Dupé 2013 for
; W3 ILC map, W7 ILC map, W9 IlC map.
;---------------------------------
; Investigates Low Quadrupole
  if keyword_set(norm) then specin = spec *1d6; ; units = muK^2
;  if not keyword_set(TheoryQuad) then TheoryQuad = 1252.d0 ; which is
;  WMAP3 theory p32
  if not keyword_set(TheoryQuad) then TheoryQuad = 1161.3421d ; for W9 data
  theory = TheoryQuad
  DataQuad = 0.
  ProbLowQuadrupole = 0.
  
  df = 5                        ; number of degrees of freedom for quadrupole = 2*l+1
  delta = theory;/2.d0/(2.d0+1)*2.d0*!dpi
  DataQuad = specin[2]
  fac = 5.d0
                                ;prob =  imsl_chisqcdf(quad*fac/delta,df,0.d0,/double)
  ProbLowQuadrupole = chisqr_pdf(DataQuad*fac/delta, df)*100.d
  Res_LowQuadrupole= {TheoryQuad: TheoryQuad,   DataQuad: DataQuad,  Prob: ProbLowQuadrupole, SimuProb:0.}
return, Res_LowQuadrupole
end
;================================================================================

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


pro s_parity, RotAlm,   lmax,  Spec, Smap, Sp, Sm
; reommanded nsiderot=64 and  nsideMap=8
 
   Alm2 = double(RotAlm.ALMre)^2 + double(RotAlm.ALMim)^2 
   NRot = RotAlm.NRot
   S = 0d
   Stest = dblarr(NRot)
;   print, Nrot
   for l= 2, Lmax do  for m=-l,l do  Stest = Stest +   (-1d)^(l+m) * reform(Alm2[l,abs(m), *] ) / Spec[l]
;   help, Stest
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

function stat_anomalies, Data,fidstat=fidstat, norm=norm, RotAlm=RotAlm, MapName=MapName, noSmap=noSmap
; Written by Anais Rassat & J.L Starck Oct 2012. 
; Modified in May 2013 by AnR & JLS: debugging and testing
;-------------------------------------------------------------------
; Details of Statistical Isotropy Tests: 
;-------------------------------------------------------------------
if not keyword_set(MapName) then MapName = 'Map'
  
 Res_LowQuad=0
 Res_QuadOctAlignment=0
 Res_PlanarOct=0
 Res_AOE=0
 Res_MirrorParity=0
                                
; ..................................................................
AnoQuadThe=0
; AnoQuadThe: set to 1 to investigate low quadrupole power
; ..................................................................
; Investigate Quadrupole/Octopole Alignment
;       Maximises equation (19) in Rassat, Starck & Dupe 2013,
;       statistic was originally defined in from de Oliveira-Costa et
;       al. 2004.
;       This version has been tested to reproduce Table 4 in Rassat, Starck & Dupe 2013, as
;       well as the favoured directions for l=2 and l-3 (see text just
;       before  section 5.3 in Rassat, Starck & Dupe 2013). Agrees
;       with Fortran code to within for nsidemap = 512 and
;       nsiderot=512.
; ..................................................................
; Investigate Planarity of Octopole
;      Maximises Equation (21) in Rassat, Starck & Dupe 2013
;      Statistic originally defined in de Oliveira-Costa et al 2004.
;      This version has been testest in May 2013 to reproduce the 't'
;      values reported in Table 5 of Rassat, Starck & Dupe 2013, which
;      were calculated using the fortran version of Stat_Anomalies. 
;      The IDL/C++ version of Stat_Anomalies agrees to nside to 3rd
;      decimal point.
;     !!! Probabilities have only been tested for nside=8 and
;     nsimu=100 !!!
; ..................................................................
;  Investigate Axis of Evil statistic
;     Results were checked against Rassat & Starck 2013 results as
;     well as against Land & Mauiejo 2006 ("AoE revisited"). The
;     results agree for Theta values with Rassat & Starck 2013, but
;     there are some discrepencies with the directions of the
;     preferred axes, especially for l=5. Results from c++/IDL code
;     seems to agree better with Land & Maguiejo 2006 (who may have
;     used the c++ version of healpix?).
; ..................................................................
;  Investigate Mirror Parity Statistic
;    Results are checked against Rassat & Starck 2013 results, Table
;    D.1
;    Probabilities (in parentheses in Table D.1 of Rassat & Starck 13
;    have not been checked yet with C++/IDL code)
; ..................................................................
;=======================   Read in map & set variables ======================
map = Data ; Map of CMB or other 
N = N_elements(map) ; Number of pixels in map
if not keyword_set(fidstat) then fidstat = set_fid_stat()
lmax = fidstat.lmaxgen
nsiderot = fidstat.nsiderot
nRot =  nsideRot^2*12L
Spec = mrs_powspec(map, lmax=20) ; calculate the Power spectrum

;==========================================================================================
;==========================================================================================
;=======================   Start Investigations of anomalies  ======================
;==========================================================================================
;==========================================================================================
if fidstat.verb eq 1 then begin
      print, '.........'
      print, 'Results for map: ', MapName
      print, 'Stat_Anomalies: Calculated rotations up to lmax = ', lmax, ' for nside = ', nsiderot
endif
;==========================================================================================
;=======================   The low quadrupole ======================
;==========================================================================================
if fidstat.lowquad.test eq 1 then begin
   Res_LowQuad= stat_lowquad(spec,theoryquad=fidstat.lowquad.theoryquad,/norm)
   if fidstat.verb eq 1 then print, 'LowQuad:', ' C(l=2)=', Res_LowQuad.DataQuad, ' (muK^2), Prob(DataQuad | TheoQuad) = ', Res_LowQuad.Prob, '(%).'
endif else Res_LowQuad = -1

; goto, DONE

;=======================   Perform Rotations  ======================
; Here we calculate the alms for nRot rotations
if fidstat.quadoct.test ne 0 and fidstat.planaroct.test ne 0 and fidstat.aoe.test ne 0 then begin
   if not keyword_set(RotAlm) then  mrs_almrot, map, RotAlm, lmax=lmax, nsideRot=nsideRot
   Alm2 = double(RotAlm.ALMre)^2 + double(RotAlm.ALMim)^2 ; calculate |a_lm(n)|^2
endif else begin
   rotalm = -1
   alm2 = -1
endelse


;==========================================================================================
;=======================  The quadrupole-octopole alignment ======================
;==========================================================================================
; Find the preferred axis for l=2 and l=3
;  to speed we can put nsideRot < 512, nsideRot=512 gives a
;precision of < 1' for the preferred axis (check in Oliveira paper)
if fidstat.quadoct.test eq 1 then begin 
; Check that the nsiderot alm for Quad/Oct is the same for the general
   if fidstat.quadoct.nsiderot ne fidstat.nsiderot then begin
      mrs_almrot, map, RotAlmQuadOct, lmax=3, nsideRot=fidstat.quadoct.nsiderot 
      Alm2QuadOct = double(RotAlm.ALMre)^2 + double(RotAlm.ALMim)^2  ; calculate |a_lm(n)|^2
   endif else begin
      RotAlmQuadOct = RotAlm
      Alm2QuadOct = Alm2
   endelse
   Res_QuadOctAlignment = stat_l2l3(map,Alm2QuadOct,RotAlmQuadOct) ; calculate stats for Quad/Oct alignment
   if fidstat.verb eq 1 then print,'Quad/Oct: ',    ' n_2 . n_3 = ', Res_QuadOctAlignment.DotProduct,  ', Angle (deg) = ', Res_QuadOctAlignment.Angle, ', Probability (%) = ', Res_QuadOctAlignment.Proba
endif else Res_QuadOctAlignment = -1
;==========================================================================================
;=======================  The planar octopole ======================
;==========================================================================================
if fidstat.planaroct.test eq 1 then begin 
   if fidstat.planaroct.nsiderot ne fidstat.nsiderot then begin
      mrs_almrot, map, RotAlmPlanarOct, lmax=3, nsideRot=fidstat.planaroct.nsiderot 
      Alm2PlanarOct = double(RotAlmPlanarOct.ALMre)^2 + double(RotAlmPlanarOct.ALMim)^2  ; calculate |a_lm(n)|^2
   endif else begin
      RotAlmPlanarOct = RotAlm
      Alm2PlanarOct = Alm2
   endelse
   Res_PlanarOct =  stat_planar_oct(spec, Alm2PlanarOct, RotAlmPlanarOct, nsimu=fidstat.planaroct.nsimu,noprob=fidstat.planaroct.noprob)
   if fidstat.verb eq 1 then print,  'Planar Oct: ', ' Stat_t = ', Res_PlanarOct.Stat_t,  ',  Probability (%) = ', Res_PlanarOct.Proba
endif else Res_PlanarOct = -1
;=======================  Evil Axis  ======================
if  fidstat.AOE.test eq 1 then begin 
   if fidstat.AOE.nsiderot ne nrot then begin
      mrs_almrot, map, RotAlmAOE, lmax=fidstat.AOE.lmax, nsideRot=fidstat.AOE.nsiderot 
      Alm2AOE = double(RotAlmAOE.ALMre)^2 + double(RotAlmAOE.ALMim)^2 ; calculate |a_lm(n)|^2
   endif else begin
      RotAlmAOE = RotAlm
      Alm2AOE = Alm2
   endelse
   Res_AOE = stat_aoe(spec, Alm2AOE, RotAlmAOE, fidstat.AOE.lmax, nsimu=fidstat.AOE.nsimu,  noprob=fidstat.AOE.noprob)
   if fidstat.verb eq 1 then print,  'Axis of Evil: ', 'Angle= ', Res_AOE.MeanAngle, ' +/- ', Res_AOE.Std, ', Proba = ', Res_AOE.Proba  
endif else Res_AOE = - 1
;=======================  Mirror Parity  ======================
if fidstat.mparity.test eq 1 then begin 
   Res_MirrorParity = stat_parity(map,NsideInMapSParity=fidstat.mparity.nside_degrade, NsideRotSParity=fidstat.mparity.nsiderot,nsimu =fidstat.mparity.nsimu,lmax=fidstat.mparity.lmax,noprob=fidstat.mparity.noprob, noSmap= noSmap)
   if fidstat.verb eq 1 then begin
      print,  'MParity: S+ = ', Res_MirrorParity.Sp,  ', Proba = ', Res_MirrorParity.ProbaSp
      print,  'Mparity: S- = ', Res_MirrorParity.Sm,  ', Proba = ', Res_mirrorParity.ProbaSm
   endif
endif else Res_MirrorParity=-1
;=======================  All Stats  ======================
DONE: 

Anomalies = { MapName: MapName, $
	                 lmax:lmax, Nsiderot:Nsiderot,  $
                     Res_LowQuad: Res_LowQuad, $  ; Quadrupole
                     Res_QuadOctAlignment: Res_QuadOctAlignment, $  ; Quad-Oct alignement
                     Res_PlanarOct: Res_PlanarOct, $         ; Oct Planarity
                     Res_AOE: Res_AOE, $         ; Axis of evil
                     Res_MirrorParity: Res_MirrorParity, $  ; Mirror Parity 
                     fidstat: fidstat $
                   }

return, Anomalies
end
