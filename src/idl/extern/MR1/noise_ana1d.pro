

PRO NOISE_ANA1D, Data, TabVar=TabVar, TabMean=TabMean, Plot=Plot, BlockSize=BlockSize, ps=ps, a0=a0, b0=b0, log=log, novar=novar,tit=tit, j=j
;+ 
; NAME: 
;       NOISE_ANA1D
;
; PURPOSE: 
;      
;       This routine separates the signal in blocks of size BlockSizexBlockSize
;       and computes inside each block the mean and the variance. The plot
;       variance versus mean gives information about the kind of noise inside
;       the data. 
;
; CALLING SEQUENCE: 
;   NOISE_ANA1D, Signal, TabVar=TabVar, TabMean=TabMean, 
;                Plot=Plot, BlockSize=BlockSize, a0=a0, b0=b0
;
; INPUTS: 
;   Image
;
; KEYED INPUT: 
;   Plot -- scalar: if set, the plot is displayed
;   BlockSize -- scalar: size of the blocks (default is 8)

;
; KEYED OUTPUT: 
;   TabVar -- IDL 1D array: array of variance
;   TabMean -- IDL 1D array: array of mean
;
; OUTPUTS: 
;   none
; 
;-
 
;------------------------------------------------------------
; parameters check
;------------------------------------------------------------
 
 IF N_PARAMS() LT 1 THEN BEGIN
   PRINT, 'CALLING SEQUENCE: ', $ 
    'NOISE_ANA1D, Signal, TabVar=TabVar, TabMean=TabMean, Plot=Plot, BlockSize=BlockSize'
   GOTO, CLOSING
 ENDIF

Image = Data
vsize = size(Image)
if vsize[0] NE 1 then BEGIN
   print, "Error: bad first parameter ..."
   print, "       Image has to be a 1 dimensional array"
   GOTO, CLOSING
 ENDIF
Nx = vsize[1]
 
if keyword_set(j) then begin 
  optmr= '-n' + STRCOMPRESS(STRING(j+1), /REMOVE_ALL)
  mr1d_trans, Image, w
  Image = w.coef(*,j-1)
end

if keyword_set(BlockSize) then SizeBlock = BlockSize else SizeBlock = 32
if not keyword_set(tit) then tit = ''
Nxb = Nx / SizeBlock
 
;------------------------------------------------------------
; function body
;------------------------------------------------------------
if keyword_set(plot) then clear

TMean = fltarr(Nxb)
TVar = fltarr(Nxb)

ind = 0L
for i=0L,Nxb-1 do BEGIN
   Debx = i*SizeBlock
   Endx = (i+1)*SizeBlock-1
   TMean[ind] = (mean(Image[Debx:Endx]))
   TVar[ind] = (sigma(Image[Debx:Endx]))
   ; TVar[ind] = get_noise(Image[Debx:Endx, Deby:Endy])
   if not keyword_set(novar) then TVar[ind] = TVar[ind]*TVar[ind]
   ind = ind + 1
   END
   
xtit='TabMean'
ytit='TabVar'
if keyword_set(log) then begin
    TMean = alog(TMean)
    TVar = alog(TVar)
    xtit='Log TabMean'
    ytit='Log TabVar'
end

ind = sort(TMean)
TabVar = TVar[ind]
TabMean = TMean[ind]

FileName='fig_'+tit+'.ps'

if keyword_set(ps) then setps, filename=FileName, /portrait

if keyword_set(plot) then BEGIN
   ;window, 1
   plot, TabMean, TabVar, xtitle=xtit, ytitle=ytit, /ynozero, psym=3,title=tit
   pause, 4
   indsort = sort(TabMean)
   n10p100 = Nxb*10/100
   Tm = (TabMean[indsort])[0:n10p100-1]
   Tv = (TabVar[indsort])[0:n10p100-1]
   ;window, 2
   ;plot, Tm, Tv, xtitle='TabMean', ytitle='TabVar', psym=3, /ynozero
   END
   
res = LINFIT( TabMean, TabVar)
a = res(1)
b =  res(0)
print, 'LINFIT(Var=a*Mean+b): a = ', res(1),  ' b = ', res(0)
legend = "a=" + strcompress(string(a,'$(f6.2)'),/remove_all) +  $
         ",b=" + strcompress(string(b,'$(f6.2)'),/remove_all)
if keyword_set(plot) then xyouts,0.6,0.8,legend,/normal
	

xx = findgen(50)/50.*(max(TabMean)-min(TabMean)) + min(TabMean) + 0.01
yy = a*xx + b
if keyword_set(plot) then begin
  ; wset, 1
  oplot,xx,yy
end
if keyword_set(ps) then endps

;print, 'V(block) = gain^2*Mean(block) + sigma^2'
;sixlin, TabMean, TabVar, a, siga, b, sigb

;print, 'gain^2 = '
;print,  b
;print,' '
;print, 'sig(gain^2) = '
;print,  sigb
;print,' '
;print, 'sigma^2 = '
;print, a
;print,' '
;print, 'sig(sigma^2) = '
;print , siga

;gain = sqrt(mean(b))
;print,' '
;print, 'gain = ', gain
;sig = sqrt(mean(a > 0))
;print, 'sigma = ', sig

;res = gain*TabMean + sig
;Err_Res = TabVar - res
;Sig_Err = sigma(Err_Res)
;index = where ( abs(Err_Res) LT Sig_Err, count)

;if (count GT 0) then BEGIN
;    TabMean = TabMean(index)
;    TabVar = TabVar(index)
;    END

;------------------------------------------------------------
; closing
;------------------------------------------------------------
 
 CLOSING:
 
  RETURN
 
 END


