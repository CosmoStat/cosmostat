;+
; NAME:
;        mrs_owtstat
;
; PURPOSE:
;	 Return statistical information relative to the bi-orthogonal wavelet transform of a given data set.
;        The return value is a 2D IDL array of 9 elements x (Number of scales-1)*(d+1) with d=0, 1, 2 for the 3 directions
;        or only d=0 if the keyword isotropic is set. The coarserst scale is not used.
;        For each scale j, we have:
;			Tab[0,j*(d+1)] = standard deviation		d = 0 : horizontal band, d = 1 : vertical band, d = 2 : diagonal band
;			Tab[1,j*(d+1)] = skewness
;			Tab[2,j*(d+1)] = Kurtosis
;			Tab[3,j*(d+1)] = Min
;			Tab[4,j*(d+1)] = Max 
;			Tab[5,j*(d+1)] = HC
;			Tab[6,j*(d+1)] = HC^+
;			Tab[7,j*(d+1)] = Cumulant order 5
;			Tab[8,j*(d+1)] = Cumulant order 6
;
;        If TabFile is set, then the statistic is computed on a set of images. Tab[*,*,f] will be the statistic related to the file TabFile[f]
;
; CALLING:
;
;      TabStat = mrs_owtstat( Data, TabStatName=TabStatName, NbrScale=NbrScale, verb=verb, TabFile=TabFile, isotropic=isotropic, 
;								survival=survival, TabSurvStat=TabSurvStat, TabAllSurvStat=TabAllSurvStat, TabSurvNu=TabSurvNu ) 
;       
; INPUTS:
;     Data -- IDL array of healpix map: Input data to analyze
;     NbrScale -- int: Number of scales. Default is 4.
;
; INPUT KEYWORDS:
;		verb: scalar -- if set, the calculated statists are printed on the screen
;		Isotropic: -- scalar, if set, directional information is not taken into account
;		TabFile -- IDL table of string, list of file where the function read the maps to be analized, in that case, on output, Imag is the last map that had been proceed and the return value is a 3D array:
;			Tab[i,j,f]	statistic i for scale j and map TabFile[f]
;		Survival -- scalar: if set, the survival function is computed instead of the different statistics
;
; OUTPUT KEYWORDS: 
;		TabStatName -- IDL table of string: TabStatName = [ "Sigma", "Skewness", "Kurtosis", "Min", "Max", "HC1", "HC2", "CUMULANT Order 5", "CUMULANT Order 6" ]
;		TabSurvStat -- 2D float array [*,j*(d+1)] survival value at scale j
;		TabAllSurvStat -- 3D double array, use together with TabFile keyword, TabAllSurvStat[*,*,f] is TabSurvStat parameter for map TabFile[f]
;		TabSurvNu -- 2D float array [*,j*(d+1)] nu survival value at scale j
;
; EXAMPLE:
;     Compute the pyramidal wavelet transform with 5 scales
;       TabStat = mrs_owtstat( Data, NbrScale=5 )
;
;       TabFile = ['File1.fits', 'File2.fits']
;       TabStat = mrs_owtstat( NbrScale=5,TabFile=TabFile )
;
; EXTERNAL CALLS:
;
; HISTORY:
;	Written: Jean-Luc Starck and Olivier Forni, 2005
;	September, 2005 File creation
;
;-----------------------------------------------------------------------------------------------------------------------------------------------------------


function mrs_owtstat, Imag, NbrScale=NbrScale, TabStatName=TabStatName, verb=verb, survival=survival, isotropic=isotropic, $
                      TabFile=TabFile, TabSurvStat=TabSurvStat, TabAllSurvStat=TabAllSurvStat, TabSurvNu=TabSurvNu

if N_PARAMS() LT 1  and not keyword_set(TabFile) then begin 
print,'CALLING SEQUENCE: TabStat = mrs_wtstat(Data, NbrScale=NbrScale, TabStatName=TabStatName, verb=verb, TabFile=TabFile, survival=survival, TabSurvStat=TabSurvStat, TabAllSurvStat=TabAllSurvStat, TabSurvNu=TabSurvNu)'
        goto, DONE
        end
	
if not keyword_set(NbrScale) then NbrScale=4

if keyword_set(TabFile) then NbrFiles = (size(TabFile))[1]

if keyword_set(TabFile) then Imag = mrs_read(TabFile[0]) $
else NbrFiles = 1

n = randomn(seed, 32,32)
Tn = get_stat(n, TabStatName=TabStatName, verb=verb)
NStat= (size(tn))[1]
TabSurvStat = 0
NSurvStat = 0
if not keyword_set(isotropic) then  NbrBand =  3*(NbrScale-1) else NbrBand = NbrScale-1
for f=0,NbrFiles-1 do begin
  if f NE 0 then Imag = mrs_read(TabFile[f])
  
  mrs_owttrans, Imag, w, NbrScale=NbrScale

  if keyword_set(survival) then begin
     NSurvStat=1001
     TabSurvStat = fltarr(NSurvStat, NbrBand)
     TabSurvNu  = fltarr(NSurvStat, NbrBand)
  end   
  TabStat = fltarr(NStat, NbrBand) 
  
 ; TabStatName = ["Sigma", "Skewness", "Kurtosis", "Min", "Max", "HC1", "HC2", "CUMULANT Order 5", "CUMULANT Order 6"]

  BorderSize = 4
  EndX = w.nx
  EndY = w.ny
  IndBand = 0
  for j =0,NbrScale-2 do begin
    HalfX = EndX/2
    HalfY = EndY/2
    size_x1 = HalfX+BorderSize
    if size_x1 LT 0 then size_x1 = 0
    size_x2 = HalfX-1-BorderSize
    if size_x2 LT 0 then size_x2 = 0
    size_x3 = EndX-1-BorderSize
    if size_x3 LT 0 then size_x3 = 0
    size_y1 = HalfY+BorderSize
    if size_y1 LT 0 then size_y1 = 0
    size_y2 = HalfY-1-BorderSize
    if size_y2 LT 0 then size_y2 = 0
    size_y3 = EndY-1-BorderSize
    if size_y3 LT 0 then size_y3 = 0
    
    BandHorizontal = w.coef[ min(size_x1,size_x3) : max(size_x1,size_x3), min(BorderSize,size_y2) : max(BorderSize,size_y2), * ];	w.coef[ size_x1 : size_x3, BorderSize : size_y2, * ]
    
    BandVertical = w.coef[ min(BorderSize,size_x2) : max(BorderSize,size_x2), min(size_y1,size_y3) : max(size_y1,size_y3), * ];	w.coef[ BorderSize : size_x2, size_y1 : size_y3, * ]
    
    BandDiagonal = w.coef[ min(size_x1,size_x3) : max(size_x1,size_x3), min(size_y1,size_y3) : max(size_y1,size_y3), * ];	w.coef[ size_x1 : size_x3, size_y1 : size_y3, * ]
    
    if keyword_set(isotropic) then begin
       BandISo = [BandHorizontal,BandVertical,BandDiagonal]
       TabStat[*,IndBand] = get_stat(BandISo, verb=verb)
       if  keyword_set(survival) then   begin
           a = survival(Band, /norm, np=Np)
           TabSurvStat[*,IndBand] = a[*,1]
	      TabSurvNu[*,IndBand] = a[*,0]
       end
       IndBand = IndBand + 1
    end else begin
      for d=0,2 do begin
         if d EQ 0 then Band = BandHorizontal else if d EQ 1 then Band =  BandVertical else Band = BandDiagonal
         TabStat[*,IndBand] = get_stat(Band, verb=verb)
	     if keyword_set(survival) then begin 
	          a = survival(Band, /norm, np=Np)
              TabSurvStat[*,IndBand] = a[*,1]
	          TabSurvNu[*,IndBand] = a[*,0]
         end
         IndBand = IndBand + 1
      end
    end
    EndX = EndX / 2
    EndY = EndY / 2
  end
    
    if keyword_set(TabFile) then BEGIN
      if f EQ 0 then BEGIN 
         vs  = size(TabStat)
         Nx = vs[1]
         Ny = vs[2]
	 TabAllStat = dblarr(nx, ny, NbrFiles)
	 if  keyword_set(survival) then   begin
	    vs  = size(TabSurvStat)
            Nx = vs[1]
            Ny = vs[2]
	    TabAllSurvStat = dblarr(nx, ny, NbrFiles)
        end
        END
       TabAllStat[*,*,f] = TabStat
       if  keyword_set(survival) then  TabAllSurvStat[*,*,f] = TabSurvStat
    END 
end


 
DONE:
  if keyword_set(TabFile) then return, TabAllStat $
  else return, TabStat
END    

