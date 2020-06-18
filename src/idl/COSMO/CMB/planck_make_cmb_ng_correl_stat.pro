;+
; NAME:
;      planck_make_cmb_ng_correl_stat
;	
; PURPOSE:
;	 Apply a set of statistics (Skewness, Kurtosis, Cumulant of order 5 and 6, min, max, Higher Criticism)  on the PLANCK CMB Map. per wavelet scale and per latiude band
;    if Mask keyword is et , the statistics at each wavelet scales are calculated only inside the mask.
;    if /Correl  is set, then  calculate the correlation per wavelet scale and per latitude band with three template  = ['Haslam 408 MHz', 'Free Free', 'IRAS 100 µm', ' 'Velocity Integrated CO Map (Dame)']
;           File names are $ISAP/data/CMB/ForegroundTemplate/408MHz_512.fits
;                                $ISAP/data/CMB/ForegroundTemplate/Halpha_256.fits"
;                                $ISAP/data/CMB/ForegroundTemplate/IRIS_100micron_mKAntenna.fits"
;                                $ISAP/data/CMB/ForegroundTemplate/lambda_wco_dht2001

;    if /Isotrop is set,  isotropy tests are performed.
;    if /Curvelet is set, curvelet statistics are calculated (Skewness, Kurtosis, Cumulant of order 5 and 6, min, max, Higher Criticism).
;    if /DWTTrans is set, Bi-orthogonal WT per Healpix face statistics are calculated  (Skewness, Kurtosis, Cumulant of order 5 and 6, min, max, Higher Criticism).
;    if /PSStat is set, Statistics inside the mask of point sources are calculted
;                                $ISAP/data/CMB/ForegroundTemplate/planck_point_sources_mask.fits"
;    if /SZStat is set, Statistics inside the mask of SZ clusters are calculted
;                                $ISAP/data/CMB/ForegroundTemplate/ERCSC_SZ_mask_rad1x5arcmin.fits"
;		
; CALLING:
;	 AllStat = planck_make_cmb_ng_correl_stat(CMBMap,  Mask=Mask, Name=Name, WTemplate_1=WTemplate_1, WTemplate_2=WTemplate_2, WTemplate_3=WTemplate_3,  WTemplate_4=WTemplate_4, NbrLat=NbrLat, $ 
;         NbrScale=NbrScale, DWTNbrScale=DWTNbrScale, CurNbrScale=CurNbrScale, Isotrop=Isotrop, correl=correl, Curvelet=Curvelet , DWTTrans=DWTTrans, PSStat= PSStat, SZStat=SZStat)
;
; INPUTS:
;   CMBMap -- Healpix Map: CMB map to analyze
;
; OUTPUTS:
;   AllStat --  StatData = {Correl1: Correl1,   --  Correlatation IDL structure relatve to template 1  (Haslam),  see mrs_wt_correlate.pro for more details
;                                    Correl2: Correl2, ,   -- Correlatation IDL structure relatve to template 2  (Free-Free)
;                                    Correl3: Correl3,    -- Correlatation IDL structure relatve to template 2   (IRAS)
;                                    StatIWT:StatIWT,  --  Undecimateded Wavelet Transform statistics. , see mrs_allstat.pro for more details
;                                    StatDWT:StatDWT,  -- Bi-orthogonal WT per Healpix face statistics
;                                    StatCUR:StatCUR,    -- Curvelet statistics
;                                    StatPS:StatPS,         -- Statistics inside the mask of point sources
;                                    StatSZ:StatSZ,        --  Statistics inside the mask of SZ clusters/
;                                    TabNameTemplate:TabNameTemplate,   -- Name of the templates
;                                    Isotropy: Iso           -- Isotropy statistics
; }  
;             
;  INPUT KEYWORDS:
;   Mask -- Healpix map  : mask of valid pixels to be used for the statistics calculation.
;   Name -- Name: Name of the CMB map
;   WTemplate_1 -- IDL Stucture: Wavelet transform of the first template
;   WTemplate_2 --  IDL Stucture: Wavelet transform of the second template
;   WTemplate_3 --  IDL Stucture: Wavelet transform of the third template
;   WTemplate_4 --  IDL Stucture: Wavelet transform of the fourth template
;   NbrLat -- number of latitude bands. Default is 25.
;   NbrScale -- int : Number of scales (default is 10) in the undecimated wavelet transform.
;	DWTNbrScale-- int : Number of scales (default is 7) in the bi-orthogonal wavelet transform per healpix face.
;	CurNbrScale-- int : Number of scales (default is 7) in the curvelet transform.
;   Isotrop -- int: if set, isotropy statistics are calculated
;   correl -- int: if set, correlation statistics are calculated
;   Curvelet -- int: if set, curvelet statistics are calculated
;   DWTTrans -- int: if set, bi-orthogonal wavelet transform statistics are calculated
;   PSStat -- int: if set, Point source  statistics are calculated
;   SZStat -- int: if set, SZ clusters  statistics are calculated
;
; EXTERNAL CALLS:
;       mrs_wt_correlate  
;   	mrs_allstat
; EXAMPLE:
;        Compute the statistics of the undecimated walelet transform and the correlation with the templates
;        Stat = planck_make_cmb_ng_correl_stat(CMBMap, Mask=Mask, /Correl, /PSStat, /SZStat)
;         
; HISTORY:
;	Written: Jean-Luc Starck, 2011
;--------------------------------------------------------------------------------------------------------



;==================================================================

function planck_make_cmb_ng_correl_stat, Ima, Mask=Mask, Name=Name, WTemplate_1=WTemplate_1, WTemplate_2=WTemplate_2, WTemplate_3=WTemplate_3,  WTemplate_4=WTemplate_4, NbrLat=NbrLat, $ 
 NbrScale=NbrScale, DWTNbrScale=DWTNbrScale, CurNbrScale=CurNbrScale, Isotrop=Isotrop, correl=correl, Curvelet=Curvelet , DWTTrans=DWTTrans, PSStat= PSStat, SZStat=SZStat, OldStat=OldStat

TabNameTemplate = ['Haslam 408 MHz', 'H-alpha', 'IRAS 100 µm', 'Velocity Integrated CO Map (Dame)']
DirISAP = getenv("ISAP")
DirTemplate = DirISAP + '/data/CMB/ForegroundTemplate/'
T1_FN = DirTemplate + '408MHz_512.fits'
T2_FN = DirTemplate +  'Halpha_256.fits'
; T3_FN = DirTemplate + "iras_100micron_2048.fits"
T3_FN = DirTemplate +  'IRIS_100micron_mKAntenna.fits'
T4_FN = DirTemplate +  "lambda_wco_dht2001.fits"
MaskSZSourcesFN = DirTemplate +  'ERCSC_SZ_mask_rad1x5arcmin.fits'
MaskSourcesFN = DirTemplate + 'planck_point_sources_mask.fits'

if not keyword_set(NbrScale) then NbrScale=10
if not keyword_set(DWTNbrScale) then DWTNbrScale=7
if not keyword_set(CurNbrScale) then CurNbrScale=7
if not keyword_set(NbrLat) then NbrLat = 25
if not keyword_set(Name) then Name = 'CMB Map'

Correl1=0
Correl2=0
Correl3=0
Correl4=0
if keyword_set(Correl) then begin
if not keyword_set(WTemplate_1) then begin
   T1 = mrs_resize(  mrs_read(T1_FN), nside=2048)
   mrs_wttrans, T1, WTemplate_1, NbrScale=NbrScale
end
if not keyword_set(WTemplate_2) then begin
   T2 = mrs_resize(  mrs_read(T2_FN), nside=2048)
   mrs_wttrans, T2, WTemplate_2, NbrScale=NbrScale
end
if not keyword_set(WTemplate_3) then begin
   T3 = mrs_resize(  mrs_read(T3_FN), nside=2048)
   mrs_wttrans, T3, WTemplate_3, NbrScale=NbrScale
end
if not keyword_set(WTemplate_4) then begin
   co = mrs_read(T4_FN)
   co = reform ( co[*,0] )
   T4 = reform(mrs_resize(co, nside=2048))
   mrs_wttrans, T4, WTemplate_4, NbrScale=NbrScale
end

  print, '   Template  ', TabNameTemplate[0]
  Correl1 = mrs_wt_correlate(Ima, T1, WT1=WT1, WT2= WTemplate_1,  NbrScale=NbrScale, Mask=Mask,  NbrLat= NbrLat)
   
  print, '   Template  ', TabNameTemplate[1]
  Correl2 = mrs_wt_correlate(Ima, T2, WT1=WT1, WT2=WTemplate_2,  NbrScale=NbrScale, Mask=Mask,  NbrLat= NbrLat)
   
  print, '   Template  ', TabNameTemplate[2]
  Correl3 = mrs_wt_correlate(Ima, T3, WT1=WT1, WT2=WTemplate_3,  NbrScale=DWTNbrScale, Mask=Mask,  NbrLat= NbrLat)
  
  print, '   Template  ', TabNameTemplate[3]
  Correl4 = mrs_wt_correlate(Ima, T4, WT1=WT1, WT2=WTemplate_4,  NbrScale=DWTNbrScale, Mask=Mask,  NbrLat= NbrLat)
  
  save, filename='xx_correl.xdr', Correl1, Correl2, Correl3, Correl4, TabNameTemplate
end ; else restore, 'xx_correl.xdr'

  print, '   WT Stat'
  StatIWT= mrs_allstat(Ima, /verb, NbrScale2D=NbrScale, Mask=Mask, /iwt)
  
  StatDWT = 0 
  if keyword_set(DWTTrans) then begin
     print, '   DWT Stat'
     StatDWT= mrs_allstat(Ima, /verb, NbrScale2D=DWTNbrScale2D, /owt)
  end
  
  StatCUR = 0
  if keyword_set(Curvelet) then begin
     print, "   CURVELET Stat"
     StatCUR = mrs_allstat(Ima, /verb, NbrScale2D=CurNbrScale, /cur)
     save, filename='xx_stat.xdr', StatIWT, StatDWT, StatCUR  
  end
  
  StatPS = 0
  if keyword_set(PSStat) then begin
     print, "   Point Source Stat"
      PSMask = mrs_read(MaskSourcesFN)
      ind = where(PSMask NE 0, c)
      PSpix = Ima[ind]
      if C GT 0 then StatPS = get_stat(PSPix, /Verb)  
      save, filename='xx_psstat.xdr', StatPS
  end
  
  StatSZ = 0
  if keyword_set(SZStat) then begin
     print, "   SZ Source Stat"
    SZMask = mrs_read(MaskSZSourcesFN)
    ind = where(SZMask NE 0, c)
    StatSZ = 0
    SZpix = Ima[ind]
    if C GT 0 then StatSZ = get_stat(SZPix, /Verb)  
    save, filename='xx_szstat.xdr', StatSZ
  end
  
  Iso=0
  if keyword_set(Isotrop) then begin
     mrs_almtrans, Ima, Alm, /complex, /tab, /norm
     Iso = mrs_isotropy(Alm=Alm, nsigma=5, TabChi2Score=TabChi2Score, TabPVal=TabPVal, TabDet=TabDet, /verb)
     save, filename='xx_iso.xdr', Iso
  end
  
  StatData = {Name:Name, Correl1: Correl1, Correl2: Correl2, Correl3: Correl3, Correl4: Correl4, StatIWT:StatIWT, StatDWT:StatDWT, StatCUR:StatCUR, StatPS:StatPS, StatSZ:StatSZ, TabNameTemplate:TabNameTemplate, Isotropy: Iso}
  OldStat = {Correl1: Correl1, Correl2: Correl2, Correl3: Correl3, StatIWT:StatIWT, StatDWT:StatDWT, StatCUR:StatCUR, StatPS:StatPS, StatSZ:StatSZ, TabNameTemplate:TabNameTemplate, Isotropy: Iso}
  save, filename='xx_allstat.xdr', StatData, OldStat

return, StatData
end

;==================================================================


