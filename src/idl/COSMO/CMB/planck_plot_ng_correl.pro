   
  
  
;==================================================
 
pro make_plot, FNVERS=FNVERS, TabStatDico, s=s, First=First, NGMax=NGMAX, Norm=Norm, Tabs=Tabs, last=last, ps=ps, TabM=TabM, DICO=DICO, TSim=TSim, Nsigma=Nsigma, abs=abs, line=line, col=col, nomean=nomean, ylog=ylog, $
right_legend= right_legend, left_legend=left_legend, top_legend=top_legend, bottom_legend= bottom_legend

; plotting stats
if not keyword_set(FNVERS) then FNVERS='fig_'
 if not keyword_set(Nsigma) then Nsigma=1
  
; restore, /verb, "stat_iwt_mask_allsimucmb.xdr"

TABMASKALLSTAT = TabStatDico

TabStatName = ["Sigma", "Skewness", "Kurtosis", "Min", "Max", "HC1", "HC2", "CUMULANT Order 5", "CUMULANT Order 6"]
TabFNStatName = ["Sigma", "Skewness", "Kurtosis", "Min", "Max", "HC1", "HC2", "Cum_Order5", "Cum_Order6"]

if not keyword_set(Tabs) then TabS=[1,2,5,6,7,8]
vs = size(tabs)
Nstat=vs[1]
window, xsize=900, ysize=500, xpos=10, ypos=600, 1
if Nstat GT 1 then  window, xsize=900, ysize=500, xpos=950,  ypos=600, 2
if Nstat GT 2 then   window, xsize=900, ysize=500, xpos=10,  ypos=100, 3 
if Nstat GT 3 then   window, xsize=900, ysize=500, xpos=950,  ypos=100, 4
if Nstat GT 4 then    window, xsize=900, ysize=500, xpos=20,  ypos=100, 5 
if Nstat GT 5 then   window, xsize=900, ysize=500, xpos=960,  ypos=100, 6


vs = size(TabStatDico)
NbrScale=vs[2]
if vs[0] EQ 2 then Nm = 1 else Nm = vs[3]
help, TabStatDico
print, ' Nm = ', Nm
print, ' Nstat = ', Nstat

if not keyword_set(TabM) then TabM = ['CCA', 'FastMEM', 'N-GenILC', 'SMICA', 'AltICA', 'LGMCA', 'N-ILC', 'Commander', 'Sevem']
if not keyword_set(DICO) then DICO = 'Wavelet'

; restore,  "stat_iwtlat_allcmb.xdr" 

; for j=0,8 do print, mean( t[2, j, *]), sigma( t[2, j, *])

if not keyword_set(s) then s=2

if not keyword_set(col) then col= indgen(Nm) + 2
if not keyword_set(line) then line= [1,2,3,4,0,0,0,0,3,4,5]
xr = indgen(NbrScale)  + 1

; for m=0,Nm-1 do begin  print, TabM[m] &  print, TABMASKALLSTAT[s,*,m]  & end

plotsimu=0
; PLOT, indgen(9), indgen(9), /YNOZERO,  XTITLE = 'Wavelet Scale', YTITLE = 'NG level', XMARGIN=[8, 8], YMARGIN=[4, 4], XSTYLE=8, /nodata
;AXIS, XAXIS=1, XTICKS=11,  XTICKN=[2250, 1125, 562, 281, 140, 70, 35, 17, 8],  XTITLE='l', XCHARSIZE = 0.7                          
; axis, xaxis=1,  xtitle='l',  charsize=1.6, xticks=9, xtickv=[2250, 1125,562, 281, 140, 70, 35, 17, 8], xrange=[2250,8] 
Freq = [2250, 1125, 562, 281, 140, 70, 35, 17, 8]
 
 if not keyword_set(First) then First=0
 if not keyword_set(Last) then Last =NbrScale-1

tek_color
!x.style=1
!y.style=1
; print, Nstat

; Nm=Nm-1
ind=0
; loadct, 39
for z = 0,Nstat-1 do begin
 
wset, z+1
s= tabS[z]
FN= FNVERS + DICO + '_' + TabFNStatName[s] + '_Scale_'  + STRC(First+1) + 'to_' +   STRC(Last+1)  
	TabMean = fltarr(NbrScale)
	TabSig = fltarr(NbrScale)
	DatatStat = TABMASKALLSTAT
if keyword_set(Tsim) then begin
	for j=0, NbrScale-1 do TabMean[j] = mean (Tsim[s,j,*])
	for j=0, NbrScale-1 do TabSig[j] = sigma (Tsim[s,j,*])
	if keyword_set(Norm) then begin
	     if keyword_set(nomean) then   for j=0, NbrScale-1 do DatatStat[s,j,*]  =  TABMASKALLSTAT[s, j,*]   / TabSig[j]  $
	     else for j=0, NbrScale-1 do DatatStat[s,j,*]  = (TABMASKALLSTAT[s, j,*] -  TabMean[j]) / TabSig[j] 
	end
	
end

if keyword_set(ps) then setps, filename=FN+'.ps'
smax = max( DatatStat[s, First:Last,*])
smin = min(DatatStat[s, First:Last,*])
smax = min([smax, 100000])
smin = max([-100000., smin])
smax = max([smax, ABS(smin)]) 
if keyword_set(NGMAX) then smax = MIN( [smax, NGMAX] )
if keyword_set(NGMAX) then  if ABS(smin) GT NGMAX then smin = - NGMAX
if keyword_set(ABS) then smin=0
Tab = reform(DatatStat[s,First:*, 0])
; plot, xtitle='Wavelet Scale', ytitle='NG level', xr[First:*],  Tab, yrange=[-smax, smax], subtitle='CMB NG: ABS(Wavelet  '+TABSTATNAME[s] + ')', line=line[0], xticks=10, xrange=[10, 0], xstyle=8, /nodata  ; , background=0, col=255

; Tit= 'CMB NG: ABS(Wavelet  '+TABSTATNAME[s] + ')'
Tit= 'NG:  ' + DICO +  '  '+ TABSTATNAME[s]  
if keyword_set(ylog) then smin=0.001
plot, xtitle='Scale', ytitle='NG level', indgen(Last-First+2)+First+1,   indgen(Last-First+1), yrange=[smin, smax], subtitle= tit, line=line[0], xticks=Last-First+2, xrange=[Last+2, First], xstyle=8, /nodata, background=255, col=0, charsize=2, ylog=ylog   ; , background=255, col=0 

; plot, xtitle='Wavelet Scale', ytitle='NG level', xr[First:*],  Tab,  title='CMB NG: ABS(Wavelet  '+TABSTATNAME[s] + ')', line=line[0], xticks=10, xrange=[0,10], xstyle=8, /ylog   ;background=255, col=0 
IE = indgen(NbrScale)+1
if keyword_set(Tsim) and not keyword_set(Norm) then  oploterror, IE[first:last],  ABS(TabMean[first:last]),  TabSig[first:last]* Nsigma, thick=2, color=0, errcol=0, line=1, errline =1 
; axis,  xaxis=1,  xtitle='l',  charsize=1., xticks=10, xtickN=reverse([4000, 2250, 1125,562, 281, 140, 70, 35, 17, 8, 0]),  xrange=[2250,8] ; , /xlog 

INDW = [-1,        0,      1,    2,     3,     4,   5,   6,   7, 8,  -1]
TN =      [4000, 2250, 1125,562, 281, 140, 70, 35, 17, 8, 0]
TN = reverse(TN[First:Last+2])
axis,  xaxis=1,  xtitle='l',  charsize=1., xticks=Last-First+2, xtickN=TN,  xrange=[2250,4],  col=0 ; , /xlog 

for m=0,Nm-1 do begin
   D = reform( (DatatStat[s,First:Last,m]))
   if keyword_set(Abs) then D = abs(D)
   oplot,  xr[First:Last], D, color=col[m], line=line[m], thick=2
  print, TabM[m], ' ', s 
 ;  print, D
end
   TL=0
    BL=0
    LL= 0
    RL=0
    if keyword_set(left_legend) then LL=1
    if keyword_set(top_legend) then TL =1
    if keyword_set(bottom_legend) then BL =1
    if keyword_set(right_legend) then RL =1

legend, TabM[0:Nm-1], linestyle= line[0:Nm-1], colors=col[0:Nm-1],  $ 
     right_legend= RL, left_legend= LL, top_legend= TL, bottom_legend= BL, $
     textcol=0, thick=3, charsize=2  ; , spacing=1, 
if keyword_set(ps) then endps $
else write_jpeg,  FN+'.jpg', tvrd(true=1),true=1  

;endps

end

DONE:
end



;==================================================

 pro make_latplot, FNVERS=FNVERS, stat=stat, TabStatDico, First=First, NGMax=NGMAX, Norm=Norm, Tabs= Tabs, j=j, Last=Last, ABS=ABS, Nsigma=Nsigma, NoMask=NoMask, ps=ps, select=select, TabM=TabM, DICO=DICO, TSim=TSim, XtabLat=XtabLat, AllScale=AllScale, line=line, col=col, nomean=nomean, win=win, ylog=ylog, $
right_legend= right_legend, left_legend=left_legend, top_legend=top_legend, bottom_legend= bottom_legend

if not keyword_set(FNVERS) then FNVERS='fig_'

if not keyword_set(Nsigma) then Nsigma=1
if not keyword_set(win) then win=0

if not keyword_set(XtabLat) then XtabLat =  86.400000 - findgen(25)*7.2
if  keyword_set(AllScale) then Abs=1

Tab = TabStatDico

TabStatName = ["Sigma", "Skewness", "Kurtosis", "Min", "Max", "HC1", "HC2", "CUMULANT Order 5", "CUMULANT Order 6"]
TabFNStatName = ["Sigma", "Skewness", "Kurtosis", "Min", "Max", "HC1", "HC2", "Cum_Order5", "Cum_Order6"]

if not keyword_set(Tabs) then TabS=[1,2,5,6,7,8]
vs = size(tabs)
Nstat=vs[1]

vs = size(TabStatDico)
NbrScale=vs[2]
if vs[0] EQ 3 then Nm = 1 else Nm = vs[4]

; help, TabStatDico
print, "Nstat = ", Nstat, ", Nmethod = ", Nm
if not keyword_set(TabM) then TabM =['GMCA', 'SMICA',  'LGMCA', 'N-ILC', 'Commander', 'Sevem']
if not keyword_set(DICO) then DICO = 'Wavelet'
 
if keyword_set(Stat) then Tabs=[stat]
 
TabSelectM = indgen(Nm)*0+1
 
; TabSelectM = [0, 0, 1, 1, 0, 1, 1, 0, 0]
   
if not keyword_set(col) then col= indgen(Nm) + 2
if not keyword_set(line) then line= [1,2,3,4,0,0,0,0,3,4,5]

xr = XtabLat
 vs = size(xr)
NbrLat = vs[1]
plotsimu=1
Freq = [2250, 1125, 562, 281, 140, 70, 35, 17, 8]

 tek_color
 ;  loadct, 39
 
!x.style=1
!y.style=1
  jj=0
 if not keyword_set(First) then First=0
 if not keyword_set(Last) then Last =NbrScale-1
if keyword_set(j) then begin
  First=j-1
  Last=j-1
  jj=j
end
if keyword_set(AllScale) then jj=0
 

TabMean = fltarr(NbrScale, NbrLat)
TabSig = fltarr(NbrScale, NbrLat)
if keyword_set(AllScale) then begin
  First=0
  Last=0
  jj=0
  TabAllMean= fltarr(NbrLat)
  TabAllSig= fltarr(NbrLat)
end

ind=1
if not keyword_set(ps) then begin
window, xsize=900, ysize=500, xpos=10, ypos=600, 1
if Nstat GT 1 then  window, xsize=900, ysize=500, xpos=950,  ypos=600, 2
if Nstat GT 2 then   window, xsize=900, ysize=500, xpos=10,  ypos=100, 3 
if Nstat GT 3 then   window, xsize=900, ysize=500, xpos=950,  ypos=100, 4
if Nstat GT 4 then    window, xsize=900, ysize=500, xpos=20,  ypos=100, 5 
if Nstat GT 5 then   window, xsize=900, ysize=500, xpos=960,  ypos=100, 6

end

for z = 0,Nstat-1 do begin
  
;  if not keyword_set(ps) then  if keyword_set(jj) then begin wset, ind  & ind=ind+1 & end
if not keyword_set(ps) then begin wset, win+ind  & ind=ind+1 & end

for j = First, Last do begin
;  if not keyword_set(jj) then begin window, ind  & ind=ind+1 & end
    s= tabS[z]
    FN= FNVERS + DICO + '_' +  'scale' + STRC(j+1) + '_' + TabFNStatName[s] + '.jpg'
    if keyword_set(AllScale) then FN= FNVERS + DICO + '_' +  'allscale' + '_' + TabFNStatName[s] + '.jpg'
    if keyword_set(ps) then begin
         FNPS = FNVERS + DICO + '_' +  'allscale' + '_' + TabFNStatName[s] + '.ps'
         setps, filename=FNPS , /portrait
    end
    
 if keyword_set(TSIM) then   begin
   for l=0, NbrLat-1 do begin
       TabMean[j,l] = mean (Tsim[s,j,l,*])
      TabSig[j,l] = sigma (Tsim[s,j,l,*])
   end
   	if keyword_set(Norm) then begin
	      if keyword_set(nomean) then   for l=0, NbrLat-1 do Tab[s,j,l,*]  = Tab[s, j,l,*] / TabSig[j,l] $
	      else  for l=0, NbrLat-1 do Tab[s,j,l,*]  = (Tab[s, j,l,*] -  TabMean[j,l] ) / TabSig[j,l] 
	end
   if keyword_set(ABS) then   TabMean[j,*] = ABS(TabMean[j,*] )
   end 
    smax = max( Tab[s, j,*,*])
    smin = min( Tab[s,j,*,*])
    smax = min([smax, 100000])
    smin = max([-100000., smin])
    smax = max([smax, ABS(smin)]) 
    if keyword_set(NGMAX) then smax = MIN( [smax, NGMAX] )
    smin=-smax
    if keyword_set(ABS) then smin = 0
    if keyword_set(ylog) then smin=0.001
 
    TabD = reform(Tab[s,j, *,0])
    
    if keyword_set(ABS) then TabD= ABS(TabD)
     if keyword_set(AllScale) then begin TabD = total(Tab[*,First:Last, *, *], 2)   & TabD = TabD[s,*,0] & end
     
    Tit='NG:  '  + DICO + ' Scale '+ STRC(j+1)+ ', l =  ' + STRC(Freq[j]) + ', ' +  TABSTATNAME[s]  
    if keyword_set(AllScale) then   Tit='NG: ABS(' + DICO + '  AllScale '+  ', l =  ' + STRC(Freq[j]) + ', ' +  TABSTATNAME[s] + ')'
        if keyword_set(ABS) then Tit='NG: ABS(' + DICO + '  Scale '+ STRC(j+1)+ ', l =  ' + STRC(Freq[j]) + ', ' +  TABSTATNAME[s] + ')'
     plot, xtitle='Latitude (degrees)', ytitle='NG level', xr,  TabD, xr=[-100,100], xticks=20, yrange=[smin, smax], subtitle=Tit, line=line[0],    /nodata, background=255, col=0 , charsize=1.5, ylog=ylog
; plot, xtitle='Wavelet Scale', ytitle='NG level', xr[First:*],  Tab,  title='CMB NG: ABS(Wavelet  '+TABSTATNAME[s] + ')', line=line[0], xticks=10, xrange=[0,10], xstyle=8, /ylog   ;background=255, col=0 
    if keyword_set(TSIM) and not keyword_set(Norm) then begin
       if keyword_set(ABS) then oploterror, xr,   ABS( reform(TabMean[j,*])),  reform( TabSig[j,*])* Nsigma, thick=2,  color=0, errcol=0 $
       else oploterror, xr,    reform(TabMean[j,*]),   reform(TabSig[j,*])*Nsigma, thick=3, color=0, errcol=0, line=1
    end  
    for m=0,Nm-1 do begin
    print, m, col[m], line[m]
        if  not keyword_set(select) or TabSelectM[m] EQ  1 then begin
          TabD = reform(Tab[s,j, *,m])
;	  info, tabD
          if keyword_set(ABS) then TabD= ABS(TabD)
   ;       if keyword_set(AllScale) then begin TabD = total( ABS(Tab[*,First:Last, *, *]), 2 )  & TabD = TabD[s,*,m] & end $
            
        if keyword_set(ABS) then oplot,  xr, TabD, color=col[m], line=line[m], thick=2 $
        else oplot,  xr, reform(  (Tab[s,j,*,m])), color=col[m], line=line[m], thick=2
        end
    end
    indm = where (TabSelectM EQ 1)
    TL=0
    BL=0
    LL= 0
    RL=0
    if keyword_set(left_legend) then LL=1
    if keyword_set(top_legend) then TL =1
    if keyword_set(bottom_legend) then BL =1
    if keyword_set(right_legend) then RL =1

    if  not keyword_set(select)  then legend, TabM[0:Nm-1], linestyle= line[0:Nm-1], colors=col[0:Nm-1],  textcol=0, charsize=2, spacing=1, thick=3, $
right_legend= RL, left_legend= LL, top_legend= TL, bottom_legend= BL  $ ;  textcol=0,  charsize=2, spacing=3
    else  legend, TabM[indm], linestyle= line[indm], colors=col[indm],  textcol=0, charsize=2, spacing=1, thick=3, $
right_legend= RL, left_legend= LL, top_legend= TL, bottom_legend= BL  ; charsize=2, , spacing=3
 
 
 ;    legend, TabM[0:Nm-1], linestyle= line[0:Nm-1], colors=col[0:Nm-1], /right_legend, /top, textcol=0, charsize=2, spacing=1, thick=3

    if keyword_set(ps) then endps $
   else write_jpeg,  FN, tvrd(true=1),true=1  

;endps
end
end

DONE:
end


;==================================================

 pro make_cumlatplot, FNVERS=FNVERS, stat=stat, TabStatDico, First=First, NGMax=NGMAX, Norm=Norm, Tabs= Tabs, j=j, Last=Last, ABS=ABS, Nsigma=Nsigma, NoMask=NoMask, ps=ps, select=select, TabM=TabM, DICO=DICO, TSim=TSim, XtabLat=XtabLat, AllScale=AllScale, line=line, col=col, nomean=nomean, win=win, ylog=ylog

if not keyword_set(FNVERS) then FNVERS='fig_'

if not keyword_set(Nsigma) then Nsigma=1
if not keyword_set(win) then win=0

if not keyword_set(XtabLat) then XtabLat =  86.400000 - findgen(25)*7.2
if  keyword_set(AllScale) then Abs=1

Tab = TabStatDico

TabStatName = ["Sigma", "Skewness", "Kurtosis", "Min", "Max", "HC1", "HC2", "CUMULANT Order 5", "CUMULANT Order 6"]
TabFNStatName = ["Sigma", "Skewness", "Kurtosis", "Min", "Max", "HC1", "HC2", "Cum_Order5", "Cum_Order6"]

if not keyword_set(Tabs) then TabS=[1,2,5,6,7,8]
vs = size(tabs)
Nstat=vs[1]

vs = size(TabStatDico)
NbrScale=vs[2]
if vs[0] EQ 3 then Nm = 1 else Nm = vs[4]

; help, TabStatDico
; print, Nstat, Nm
if not keyword_set(TabM) then TabM = ['CCA', 'FastMEM', 'N-GenILC', 'SMICA', 'AltICA', 'LGMCA', 'N-ILC', 'Commander', 'Sevem']
if not keyword_set(DICO) then DICO = 'Wavelet'
 
if keyword_set(Stat) then Tabs=[stat]
 
 TabSelectM = indgen(Nm)*0+1
 
; TabSelectM = [0, 0, 1, 1, 0, 1, 1, 0, 0]
   
if not keyword_set(col) then col= indgen(Nm) + 2
if not keyword_set(line) then line= [1,2,3,4,0,0,0,0,3,4,5]

xr = XtabLat
 vs = size(xr)
NbrLat = vs[1]
plotsimu=1
Freq = [2250, 1125, 562, 281, 140, 70, 35, 17, 8]

 tek_color
 ;  loadct, 39
 
!x.style=1
!y.style=1
  jj=0
 if not keyword_set(First) then First=0
 if not keyword_set(Last) then Last =NbrScale-1
if keyword_set(j) then begin
  First=j-1
  Last=j-1
  jj=j
end
if keyword_set(AllScale) then jj=0
 

TabMean = fltarr(NbrScale, NbrLat)
TabSig = fltarr(NbrScale, NbrLat)
MaxTab = fltarr(NbrLat, Nm)
 
if keyword_set(AllScale) then begin
  First=0
  Last=0
  jj=0
  TabAllMean= fltarr(NbrLat)
  TabAllSig= fltarr(NbrLat)
end

ind=1
if not keyword_set(ps) then begin
window, xsize=900, ysize=500, xpos=10, ypos=600, 1
if Nstat GT 1 then  window, xsize=900, ysize=500, xpos=950,  ypos=600, 2
if Nstat GT 2 then   window, xsize=900, ysize=500, xpos=10,  ypos=100, 3 
if Nstat GT 3 then   window, xsize=900, ysize=500, xpos=950,  ypos=100, 4
if Nstat GT 4 then    window, xsize=900, ysize=500, xpos=20,  ypos=100, 5 
if Nstat GT 5 then   window, xsize=900, ysize=500, xpos=960,  ypos=100, 6

end

for z = 0,Nstat-1 do begin

;  if not keyword_set(ps) then  if keyword_set(jj) then begin wset, ind  & ind=ind+1 & end
if not keyword_set(ps) then begin wset, win+ind  & ind=ind+1 & end

for j = First, Last do begin
;  if not keyword_set(jj) then begin window, ind  & ind=ind+1 & end
    s= tabS[z]
    FN= FNVERS + DICO + '_' +  'cumlat_scale' + STRC(j+1) + '_' + TabFNStatName[s] + '.jpg'
    if keyword_set(AllScale) then FN= FNVERS + DICO + '_' +  'allscale' + '_' + TabFNStatName[s] + '.jpg'
    if keyword_set(ps) then begin
         FNPS = FNVERS + DICO + '_' +  'allscale' + '_' + TabFNStatName[s] + '.ps'
         setps, filename=FNPS , /portrait
    end
    
 if keyword_set(TSIM) then   begin
   for l=0, NbrLat-1 do begin
       TabMean[j,l] = mean (Tsim[s,j,l,*])
      TabSig[j,l] = sigma (Tsim[s,j,l,*])
   end
   	if keyword_set(Norm) then begin
	      if keyword_set(nomean) then   for l=0, NbrLat-1 do Tab[s,j,l,*]  = Tab[s, j,l,*] / TabSig[j,l] $
	      else  for l=0, NbrLat-1 do Tab[s,j,l,*]  = (Tab[s, j,l,*] -  TabMean[j,l] ) / TabSig[j,l] 
	end
   if keyword_set(ABS) then   TabMean[j,*] = ABS(TabMean[j,*] )
   end 
    smax = max( Tab[s, j,*,*])
    smin = min( Tab[s,j,*,*])
    smax = min([smax, 100000])
    smin = max([-100000., smin])
    smax = max([smax, ABS(smin)]) 
    if keyword_set(NGMAX) then smax = MIN( [smax, NGMAX] )
    smin=-smax
    if keyword_set(ABS) then smin = 0
    if keyword_set(ylog) then smin=0.001

    TabD = reform(Tab[s,j, *,0])
    if keyword_set(ABS) then TabD= ABS(TabD)
     if keyword_set(AllScale) then begin TabD = total(Tab[*,First:Last, *, *], 2)   & TabD = TabD[s,*,0] & end
     
    Tit='NG:  '  + DICO + ' Scale '+ STRC(j+1)+ ', l =  ' + STRC(Freq[j]) + ', ' +  TABSTATNAME[s]  
    if keyword_set(AllScale) then   Tit='NG: ABS(' + DICO + '  AllScale '+  ', l =  ' + STRC(Freq[j]) + ', ' +  TABSTATNAME[s] + ')'
        if keyword_set(ABS) then Tit='NG: ABS(' + DICO + '  Scale '+ STRC(j+1)+ ', l =  ' + STRC(Freq[j]) + ', ' +  TABSTATNAME[s] + ')'
      plot, xtitle='Fsky (%)', ytitle='NG level', xr,  TabD[0:NbrLat/2-1], xr=[100,10], xticks=10, yrange=[smin, smax], subtitle=Tit, line=line[0],    /nodata, background=255, col=0 , charsize=1.5, ylog=ylog
  ;  plot, xtitle='ABS(Minimum Latitude) (degrees)', ytitle='NG level', xr[NbrLat/2:*],  TabD[NbrLat/2:*], xr=[100,0], xticks=20, yrange=[smin, smax], subtitle=Tit, line=line[0],    /nodata, background=255, col=0 , charsize=1.5, ylog=ylog
; plot, xtitle='Wavelet Scale', ytitle='NG level', xr[First:*],  Tab,  title='CMB NG: ABS(Wavelet  '+TABSTATNAME[s] + ')', line=line[0], xticks=10, xrange=[0,10], xstyle=8, /ylog   ;background=255, col=0 
    if keyword_set(TSIM) and not keyword_set(Norm) then begin
       if keyword_set(ABS) then oploterror, xr,   ABS( reform(TabMean[j,*])),  reform( TabSig[j,*])* Nsigma, thick=2,  color=0, errcol=0 $
       else oploterror, xr,    reform(TabMean[j,*]),   reform(TabSig[j,*])*Nsigma, thick=3, color=0, errcol=0, line=1
    end  
    for m=0,Nm-1 do begin
        if  not keyword_set(select) or TabSelectM[m] EQ  1 then begin
          TabD = reform(Tab[s,j, *,m])
          if keyword_set(ABS) then TabD= ABS(TabD)
   ;       if keyword_set(AllScale) then begin TabD = total( ABS(Tab[*,First:Last, *, *]), 2 )  & TabD = TabD[s,*,m] & end $
            
        if keyword_set(ABS) then oplot, xr,  TabD[0:NbrLat/2-1], color=col[m], line=line[m], thick=2 $
        else oplot,  xr, reform(  (Tab[s,j,0:NbrLat/2-1,m])), color=col[m], line=line[m], thick=2
        end
    end
    indm = where (TabSelectM EQ 1)
    if  not keyword_set(select)  then legend, TabM[0:Nm-1], linestyle= line[0:Nm-1], colors=col[0:Nm-1], /right_legend, /top, textcol=0, charsize=2, spacing=1, thick=3  $ ;  textcol=0,  charsize=2, spacing=3
    else  legend, TabM[indm], linestyle= line[indm], colors=col[indm], /right_legend, /top, textcol=0, charsize=2, spacing=1, thick=3  ; charsize=2, , spacing=3
 
 
 ;    legend, TabM[0:Nm-1], linestyle= line[0:Nm-1], colors=col[0:Nm-1], /right_legend, /top, textcol=0, charsize=2, spacing=1, thick=3

    if keyword_set(ps) then endps $
   else write_jpeg,  FN, tvrd(true=1),true=1  

;endps
end  ; end j 
 
 ; TabD = reform(Tab[s,j, lat, Method])
 PlotMaxScale=0
 if keyword_set(PlotMaxScale) then begin
for l=0, NbrLat-1 do begin
T1 = reform(Tab[s,*, *,*])
 if keyword_set(ABS) then T1= ABS(T1)
  for m=0,Nm-1 do  MaxTab[l,m] = max( T1[*,l,m])
 end
 smin = min(MaxTab)
 smax= max(MaxTab)
 TabD = reform(MaxTab[*,0])
 
    FN= FNVERS + DICO + '_' +  'cumlat_maxscale' +  '_' + TabFNStatName[s] + '.jpg'
    if keyword_set(ps) then begin
         FNPS = FNVERS + DICO + '_' +  'cumlat_maxscale' + '_' + TabFNStatName[s] + '.ps'
         setps, filename=FNPS , /portrait
    end
 
 Tit='NG:  '  + DICO + ' Max over Scales ' +  TABSTATNAME[s]  
  if keyword_set(ABS) then Tit='Max of  ABS(' + TABSTATNAME[s] + ')' +  ' over ' + DICO  + ' Scales versus Mask size'
     plot, xtitle='Mask size: ABS(Minimum Latitude) (degrees)', ytitle='NG level', xr,  TabD[0:NbrLat/2-1], xr=[100,10], xticks=20, yrange=[smin, smax], subtitle=Tit, line=line[0],    /nodata, background=255, col=0 , charsize=1.5, ylog=ylog
  for m=0,Nm-1 do begin
          TabD = reform(MaxTab[*,m])
          if keyword_set(ABS) then TabD= ABS(TabD)
          oplot,  xr, TabD, color=col[m], line=line[m], thick=2 
         end
;  oplot,  col=0, xr, TabD*0+3, line=1
;  oplot,  col=0, xr, TabD*0+5, line=2

     indm = where (TabSelectM EQ 1)
    if  not keyword_set(select)  then legend, TabM[0:Nm-1], linestyle= line[0:Nm-1], colors=col[0:Nm-1], /right_legend, /top, textcol=0, charsize=2, spacing=1, thick=3  $ ;  textcol=0,  charsize=2, spacing=3
    else  legend, TabM[indm], linestyle= line[indm], colors=col[indm], /right_legend, /top, textcol=0, charsize=2, spacing=1, thick=3  ; charsize=2, , spacing=3
 
   if keyword_set(ps) then endps $
   else write_jpeg,  FN, tvrd(true=1),true=1  
end ; keyword_set(PLotMaxScale)
end  ;end stat

DONE:
end


;==================================================
 
pro make_correl_plot, FNVERS=FNVERS, TabCorrel, temp=temp, ymax=ymax, abs=abs, Norm=Norm, First=First, Last=Last, ps=ps, old=old, nsig=nsig, Method=Method, TabSim=TabSim , line=line, col=col, nomean=nomean, ylog=ylog, NGMAX=NGMAX, $
right_legend= right_legend, left_legend=left_legend, top_legend=top_legend, bottom_legend= bottom_legend


if not keyword_set(FNVERS) then FNVERS='fig_'
if not keyword_set(nsig) then nsig =1.

if not keyword_set(temp) then temp=1
t=temp-1

thisLetter = "154B
muLetter = '!4' + String(thisLetter) + '!X'
TabNameTemplate = ['Haslam 408 MHz', 'Free Free', 'IRAS 100 '+muletter + 'm', 'H1']
FNNameTemplate = ['Haslam', 'FreeFree', 'IRAS', 'H1']

vs = size(TabCorrel)
NbrScale=vs[1]
if vs[0] EQ 1 then Nm = 1 $
else Nm = vs[2]

if not keyword_set(Method) then Method = ['CCA', 'FastMEM', 'N-GenILC', 'SMICA', 'AltICA', 'LGMCA', 'N-ILC', 'Commander', 'Sevem']

vs = size(FNNameTemplate)
NbrTemp = vs[1]
 
if not keyword_set(col) then col= indgen(Nm) + 2
if not keyword_set(line) then line= [1,2,3,4,0,0,0,0,3,4,5]
xr = indgen(NbrScale)  + 1
Freq = [2250, 1125, 562, 281, 140, 70, 35, 17, 8, 4]


 if not keyword_set(First) then First=0
 if not keyword_set(Last) then Last =NbrScale-1
if keyword_set(j) then begin
  First=j
  Last=j
end

Tab = TabCorrel
 
TabM = Method
FN= FNVERS +  'correl_' + FNNameTemplate[t] + '_Scale_'  + STRC(First+1) + 'to_' +   STRC(Last+1)  
if keyword_set(ps) then setps, filename=FN+'.ps'

TabMean = fltarr(NbrScale)
TabSig = fltarr(NbrScale)
if keyword_set(ABS) then Tab = ABS(Tab)

if keyword_set(TabSIM) then begin
for j=0, NbrScale-1 do TabMean[j] = mean (TabSIM[j,*])
for j=0, NbrScale-1 do TabSig[j] = sigma (TabSIM[j,*])
	if keyword_set(Norm) then begin
     	if keyword_set(nomean) then  for j=0, NbrScale-1 do Tab[j,*]  = Tab[j,*]  / TabSig[j]  $
	    else for j=0, NbrScale-1 do Tab[j,*]  = (Tab[j,*] -  TabMean[j]) / TabSig[j] 
	end
end

 if not keyword_set(ymax) and not keyword_set(Norm) then ymax = 1 $
 else ymax = max(Tab)
 ymin=-ymax
 if   keyword_set(abs) then ymin = 0


if keyword_set(NGMAX) then NGMAXCorrel = NGMAX else NGMAXCorrel = ymax

    smax = max(Tab[First:Last,*])
    smin = min( Tab[First:Last,*])
    smax = max([smax, ABS(smin)]) 
    smax = MIN( [smax, NGMAXCorrel] )
    smin=-smax
    if keyword_set(ABS) then smin = 0

print, first, last
Tit='  Correlation: Wavelet  '+ TabNameTemplate[t] 
if keyword_set(ABS) then Tit='  Correlation: ABS(Wavelet  '+ TabNameTemplate[t] + ')'

plot, xtitle='Wavelet Scale', ytitle='Correlation', indgen(Last-First+2)+First+1,   indgen(Last-First+1), yrange=[smin, smax], subtitle=Tit, line=line[0], xticks=Last-First+2, xrange=[Last+2, First], xstyle=8, /nodata, background=255, col=0 , charsize=2, ylog=ylog

; plot, xtitle='Wavelet Scale', ytitle='Correlation', xr[First:Last],  Tab[First:Last,0], yrange=[smin, smax], subtitle='  Correlation: ABS(Wavelet  '+ TabNameTemplate[t] + ')', line=line[0], xticks=Last-First+1, xrange=[Last, First], xstyle=8, /nodata   ;background=255, col=0 
; oploterror, indgen(9)+1,  ABS(TabMean),  TabSig, thicka
INDW = [-1,        0,      1,    2,     3,     4,   5,   6,   7, 8, 9, -1]
TN =      [4000, 2250, 1125,562, 281, 140, 70, 35, 17, 8, 4,0]
TN = reverse(TN[First:Last+2])
; TN = reverse(TN[First:Last])

axis,  xaxis=1,  xtitle='l',  charsize=1., xticks=Last-First+2, xtickN=TN,  xrange=[2250,4] , color=0; , /xlog 
if keyword_set(TabSIM) and not keyword_set(Norm) then oploterror, indgen(NbrScale)+1,  ABS(TabMean),  TabSig*nsig, thick=2, color=0, errcol=0, line=1
for m=0,Nm-1 do oplot,  xr[First:Last], reform( (Tab[First:Last,m])), color=col[m], line=line[m], thick=2

  TL=0
    BL=0
    LL= 0
    RL=0
    if keyword_set(left_legend) then LL=1
    if keyword_set(top_legend) then TL =1
    if keyword_set(bottom_legend) then BL =1
    if keyword_set(right_legend) then RL =1

legend, TabM[0:Nm-1], linestyle= line[0:Nm-1], colors=col[0:Nm-1],  spacing=3,  textcol=0, thick=3, charsize=2, $
right_legend= RL, left_legend= LL, top_legend= TL, bottom_legend= BL  ; , spacing=1, 

if keyword_set(ps) then endps $
else write_jpeg,  FN + '.jpg', tvrd(true=1),true=1  

end

;==================================================

pro make_lat_correl_plot, FNVERS=FNVERS,  TabCorrel,  ymax=ymax, abs=abs,  Norm=Norm, First=First, Last=Last, j=j,  Nsigma=Nsigma, ps=ps, temp=temp, Method=Method, XtabLat=XtabLat, TabSim=TabSim, line=line, col=col, nomean=nomean, win=win


if not keyword_set(FNVERS) then FNVERS = 'fig_'
if not keyword_set(XtabLat) then XtabLat =  86.400000 - findgen(25)*7.2
if not keyword_set(Nsigma) then Nsigma=1
if not keyword_set(win) then win=0

if not keyword_set(t) then t=0
TabNameTemplate = ['Haslam 408 MHz', 'Free Free', 'IRAS 100 µm', 'H1']
FNNameTemplate = ['Haslam', 'FreeFree', 'IRAS', 'H1']

if not keyword_set(temp) then temp=1
t=temp-1

vs = size(TabCorrel)
NbrScale=vs[1]
print, vs
help, Nm
if vs[0] EQ 2 then Nm = 1 else Nm = vs[3]

if not keyword_set(col) then col= indgen(Nm) + 2
if not keyword_set(line) then line= [1,2,3,4,0,0,0,0,3,4,5]
; col[8] = 4

xr = XtabLat
vs = size(xr)
NbrLat = vs[1]

Freq = [2250, 1125, 562, 281, 140, 70, 35, 17, 8, 4]

 if not keyword_set(First) then First=0
 if not keyword_set(Last) then Last =NbrScale-1
if keyword_set(j) then begin
   First=j-1
  Last=j-1
  jj=j
end
print, first, last

;;TabMean = fltarr(NbrScale, NbrLat)
;TabSig = fltarr(NbrScale, NbrLat)

 tek_color
!x.style=1
!y.style=1

Tab =  TabCorrel
TabM = Method
 
FN= FNVERS + 'correl_' + FNNameTemplate[t] + '.jpg'

TabMean = fltarr(NbrScale, NbrLat)
TabSig = fltarr(NbrScale, NbrLat)

if keyword_set(ABS) then Tab = ABS(Tab)
;print, "TABM", TabM
;print, col

ind=0
j=0
for j = First, Last do begin
; print, j
    window, ind+ win
    ind=ind+1
    if keyword_set(TabSIM) then begin
    for l=0, NbrLat-1 do  TabMean[j,l] = mean (TabSIM[j, l, *])
    for l=0, NbrLat-1 do  TabSig[j, l] = sigma (TabSIM[j,l, *])
    if keyword_set(Norm) then begin
	     if keyword_set(nomean) then for l=0, NbrLat-1 do Tab[j,l,*]  = Tab[j,l,*]  / TabSig[j,l] $
	     else for l=0, NbrLat-1 do Tab[j,l,*]  = (Tab[j,l,*] -  TabMean[j,l] ) / TabSig[j,l] 
	end
    end
    help, ymax
     if not keyword_set(ymax)    then  ymax = max(Tab)
    ymin=-ymax
    if   keyword_set(abs) then ymin = 0


    FN = FNVERS + 'scale' + STRC(j+1) + '_' + FNNameTemplate[t]  

    if keyword_set(ps) then begin
          FNPS =FN + '.ps'
          setps, filename=FNPS , /portrait
    end
    NGMAX= ymax
    smax = max(Tab[j,*,*])
    smin = min( Tab[j,*,*])
    smax = max([smax, ABS(smin)]) 
    smax = MIN( [smax, NGMAX] )
    smin=-smax
    if keyword_set(ABS) then smin = 0
       
     Tit='Correlation of Wavelet  Scale '+ STRC(j+1)+ ', l =  ' + STRC(Freq[j]) + ' with ' +  FNNameTemplate[t] + ')'
     plot, xtitle='Latitude (degrees)', ytitle='Correlation', xr,  Tab[j,*,0], xr=[-100,100], xticks=20, yrange=[smin, smax], subtitle=Tit, line=line[0],    /nodata , background=255, col=0, thick=2, charsize=2   ;background=255, col=0 

     me = reform(TabMean[j,*])
     se = reform (TabSig[j,*])
     if keyword_set(TabSIM) and not keyword_set(Norm) then begin
        if keyword_set(ABS) then oploterror, xr,   ABS( me),  se * Nsigma, thick=2, color=0, errcol=0, line=1, errline =1 $
        else oploterror, xr,    me,   se *Nsigma, thick=3, color=0, errcol=0, line=1, errline=1
    end
    
  ;  help, Tab
    TabD =  Tab
    if keyword_set(ABS) then TabD = ABS(TabD)
    for m=0,Nm-1 do  oplot,  xr, reform( TabD[j,*,m]), color=col[m], line=line[m], thick=3
    legend, TabM[0:Nm-1], linestyle= line[0:Nm-1], colors=col[0:Nm-1], /right_legend, /top, textcol=0, thick=3, charsize=2.5, spacing=1
    if keyword_set(ps) then endps $
    else begin
     print, "JPEG ", FN+ '.jpg'
     write_jpeg,  FN+ '.jpg', tvrd(true=1),true=1  
     end
end

DONE:

end

;==================================================

pro make_iso_plot,  TabIso, xr=xr, SimTab=SimTab, Method=Method, line=line, col=col

PRE = 'ISO/iso_inp_masked_cmb_'

SIM=0
if keyword_set(SIM) then begin
Nf=20

   if  keyword_set(SimTab) then  SimTabISO = SimTab
   for f=0,Nf-1 do begin
     FN = PRE + strc(f+1)+'.xdr'
     if not keyword_set(SimTab) then begin
        restore, /verb, FN
        if f eq 0 then SimTabISO = replicate(s, Nf)
        SimTabISO[f] = s
      end
  end
if  not keyword_set(SimTab) then  begin
    SimTab = SimTabISO
    save, filename='iso_allcmbxdr', SimTab
end
NbrL=3000
TS = SimTab.TABPVAL
TabMean = fltarr(NbrL)
TabSig = fltarr(NbrL)
 for l = 0, NbrL-1 do begin
      TabMean[l] = mean (TS[l, *])
      TabSig[l] = sigma (TS[l, *])
end
end

TabM = Method
vs = size(TabM)
Nf = vs[1]

if not keyword_set(xr) then xr=[0,3000]
for f=0,Nf-1 do begin
    window, f
    print, f, ' ', TabM[f], TabISO[f].pval
    plot, TabISO[f].TABPVAL, /ynoz, tit=TabM[f], yrange=[0.5, 1], xrange=xr
    if keyword_set(SIM) then  begin
         oplot, TabMean-TabSig, line=1
         oplot, TabMean-TabSig*2, line=2
    end
end
end


;==================================================

 
