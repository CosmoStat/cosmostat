function cmb_anomalies, Data, lmax=lmax, nside=nside, TheoryQuadrupole=TheoryQuadrupole, norm=norm, verb=verb, AnoQuadThe=AnoQuadThe,  AnoQuadOctAlignment=AnoQuadOctAlignment,  AnoPlanarOctopole=AnoPlanarOctopole, RotAlm=RotAlm
; Written by Anais Rassat & J.L Starck Oct 2012. 
;  
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
N = N_elements(map)

if not keyword_set(lmax) then lmax = 5L
if keyword_set(Verb) then print, 'Lmax = ',  lmax
if keyword_set(norm) then map = map / 1d6/2.725
if not keyword_set(TheoryQuadrupole) then TheoryQuadrupole = 1252.d0 ; which is WMAP3 theory p32
if not keyword_set(nside) then nside = 32

 Spec = mrs_powspec(map, lmax=20)

;=======================   The low quadrupole ======================
AnoQuadThe=0

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
end

if not keyword_set(RotAlm) then  mrs_almrot, map, RotAlm, lmax=lmax
 Alm2 = RotAlm.ALMre^2 + RotAlm.ALMim^2 

;=======================  The quadrupole-octopole alignment ======================

AnoQuadOctAlignment = 1

TabQl=1
; Find the preferred axis for l=2 and l=3
;  to speed we can put nside < 512, nside=512 gives a precision of < 1' for the preferred axis (check in Oliveira paper)
if  keyword_set(AnoQuadOctAlignment) then begin 
    TabSum = dblarr(N)
    TabQl = dblarr(lmax+1, 5)
    for l=2,lmax do begin
      TabL = reform(Alm2[l,0:l, *])
      indm = lindgen(l+1)
      indm2 = indm * indm
      for i=0L,N-1 do TabSum[i] = 2. * total( TabL[*,i] * indm2 )
      TabQl[l,0] = max(TabSum, indAngle)
      TabQl[l,1]  = RotAlm.THETAPHI[indAngle, 0]
      TabQl[l,2]  = RotAlm.THETAPHI[indAngle, 1]
      ang2lb, TabQl[l,1],  TabQl[l,2] , Coord_l,  Coord_b
      TabQl[l,3]  = Coord_l
      TabQl[l,4]  = Coord_b
   end
   if keyword_set(Verb) then for l=2,lmax do  print, 'l= ', strc(l), ' Ql = ', strc(TabQl[l,0]), ', Angle (l,b) = (', strc(TabQl[l,3]) , ', ', strc(TabQl[l,4]) , ')'
end

;=======================  The planar octopole ======================
AnoPlanarOctopole=1

M3Coord_l=0
M3Coord_b=0
StatM3=0
if  keyword_set(AnoPlanarOctopole) then begin 
    l=3
    TabL = reform(Alm2[l,0:l, *])
    TabM3 = dblarr(N)
    for i=1L,N-1 do TabM3[i] =  ( 2. * TabL[3,i] )  /  ( 2. * total( TabL[1:*,i]  ) + TabL[0,i]  )
    StatM3 = max(TabM3, indAngle)
    IndM3 = indAngle
    M3Theta  = RotAlm.THETAPHI[IndM3, 0]
    M3Phi  = RotAlm.THETAPHI[IndM3, 1]
    ang2lb, TabQl[l,1],  TabQl[l,2] , M3Coord_l,  M3Coord_b
 end

;=======================  Evil Axis  ======================

AnoAxisOfEvil=0

TabRmaxMN= 0
 TabMaxIndAngle= 0
 TabMaxIndM= 0
 TabRmaxTheta= 0
 TabRmaxPhi= 0
 TabRmaxLCoord= 0
 TabRmaxRCoordi = 0
 if  keyword_set(AnoAxisOfEvil) then begin 
    ll = lindgen(lmax+1)
    ll = (2*ll+1)
    ClDen =  ll * Spec
    TabMax = dblarr(N)
    TabIndM= dblarr(N)
    TabRmaxMN = dblarr(lmax+1)
    TabMaxIndM = dblarr(lmax+1)
    TabMaxIndAngle = dblarr(lmax+1)
    TabRmaxTheta = dblarr(lmax+1)
    TabRmaxPhi = dblarr(lmax+1)
    TabRmaxLCoord = dblarr(lmax+1)
    TabRmaxRCoordi = dblarr(lmax+1)

    for l=2,lmax do begin
      TabL = reform(Alm2[l,0:l, *])
      Clm = TabL
      Clm[1:*, *] = 2. * Clm[1:*, *]
      Clm = Clm / ClDen[l]
      for i=0L,N-1 do begin
          TabMax[i] = max( Clm[*,i], indm )
          TabIndM[i] = indm
      end
      TabRmaxMN[l] = max( TabMax, indi )
      TabMaxIndAngle[l] = indi
      TabMaxIndM[l] = TabIndM[indi]  
      TabRmaxTheta[l]  = RotAlm.THETAPHI[indi, 0]
      TabRmaxPhi[l]  = RotAlm.THETAPHI[indi, 1]
      ang2lb, RotAlm.THETAPHI[indi, 0],  RotAlm.THETAPHI[indi, 1] , RCoord_l,  RCoord_b     
      TabRmaxLCoord[l]  = RCoord_l
      TabRmaxBCoord[l]  = RCoord_b
      end
      if keyword_set(Verb) then print, 'R_maxMN [l] = ', TabRmaxMN[2:*]
end


;=======================  The S-Stat  ======================

Anomalies = { lmax:lmax, N:N, $
                     TheoryQuadrupole: TheoryQuadrupole, $  ; Quadrupole
                     DataQuadrupole: DataQuadrupole, $
                     ProbLowQuadrupole: ProbLowQuadrupole, $  
                     QuadOctAlignment_Ql: TabQl, $   ; Quad-Oct alignement
                     StatPlanarOctopole: StatM3, $
                     M3Coord_l: M3Coord_l, $         ; Oct Planarity
                     M3Coord_b: M3Coord_b, $
                      TabRmaxMN: TabRmaxMN, $    ;  Axis of evil
                      TabMaxIndAngle: TabMaxIndAngle, $
                      TabMaxIndM: TabMaxIndM, $
                      TabRmaxTheta: TabRmaxTheta, $
                      TabRmaxPhi: TabRmaxPhi, $
                      TabRmaxLCoord: TabRmaxLCoord, $
                      TabRmaxRCoordi: TabRmaxRCoordi}


return, Anomalies
end
