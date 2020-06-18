

PRO NOISE_ANA, Data, DataNorm=DataNorm, ImaNoise=ImaNoise, Getnoise=GetNoise, TabVar=TabVar, TabMean=TabMean, Plot=Plot, BlockSize=BlockSize, ps=ps, a0=a0, b0=b0, $
               log=log, novar=novar,nofit=nofit, tit=tit, j=j, yrange=yrange, xrange=xrange
;+ 
; NAME: 
;       NOISE_ANA
;
; PURPOSE: 
;      
;       This routine separates the image in blocks of size BlockSizexBlockSize (8 by default)
;       and computes inside each block the mean and the variance. The plot
;       variance versus mean gives information about the kind of noise inside
;       the data. If /Getnoise is set, the  local variance of the data calculation is replaced
;       by the local estimated variance of the noise. 
;
; CALLING SEQUENCE: 
;
;   NOISE_ANA, Image, TabVar=TabVar, TabMean=TabMean, 
;                Plot=Plot, BlockSize=BlockSize, a0=a0, b0=b0
;
; INPUTS: 
;   Image
;
; KEYED INPUT: 
;   Plot -- scalar: if set, the plot is displayed
;   BlockSize -- scalar: size of the blocks (default is 8)
;   Getnoise -- scalar: if set, replace the local variance of the data by the local estimated variance of the noise
;   PS -- scalar: if set, plot the result in the postcript file
;   tit -- string: filename of the postcript file = 'fig_' + tit + '.ps'
;   log -- scalar: if set, plot the LOG(Variance) versus LOG(Mean).
;   novar -- scale: if set, plot the STD versus the Mean instead of the Variance versus Mean
;   j -- integer: if set, the noise analysis method is applied on the jth scale of the isotropic wavelet transform
;
;
; KEYED OUTPUT: 
;   TabVar -- IDL 1D array: array of variance
;   TabMean -- IDL 1D array: array of mean
;   ImaNoise -- IDL 2D array:  image of the local estimation of the standard deviation
;   DataNorm -- IDL 2D array:  normalized data : DataNorm = Data / ImaNoise
;   a0 -- IDL float: result of the line fit on the curve TabVar(TabMean):  
;   b0 -- IDL float:                   TabVar = a0 * TabMean + b0
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
    'NOISE_ANA, Image, TabVar=TabVar, TabMean=TabMean, j=j, Plot=Plot, log=log, novar=novar, BlockSize=BlockSize(def=8), GetNoise=GetNoise, DataNorm=DataNorm, ImaNoise=ImaNoise, nofit=nofit, tit=tit, j=j, ps=ps'
   GOTO, CLOSING
 ENDIF

Image = float(Data)
vsize = size(Image)
if vsize[0] NE 2 then BEGIN
   print, "Error: bad first parameter ..."
   print, "       Image has to be a 2 dimensional array"
   GOTO, CLOSING
 ENDIF
Nx = vsize[1]
Ny = vsize[2]

if keyword_set(j) then begin 
  optmr= '-n' + STRCOMPRESS(STRING(j+1), /REMOVE_ALL)
  mr_transform, Image, w
  Image = w(*,*,j-1)
end

if keyword_set(BlockSize) then SizeBlock = BlockSize else SizeBlock = 8
if not keyword_set(tit) then tit = ''
Nxb = Nx / SizeBlock
Nyb = Ny / SizeBlock
if Nxb*SizeBlock NE Nx then Nxb = Nxb+1
if Nyb*SizeBlock NE Ny then Nyb = Nyb+1

;------------------------------------------------------------
; function body
;------------------------------------------------------------
if keyword_set(plot) then clear

TMean = fltarr(Nxb*Nyb)
TVar = fltarr(Nxb*Nyb)
ImaSigma = fltarr(Nxb,Nyb)
DataNorm = Image
ind = 0L
for i=0,Nxb-1 do BEGIN
for j=0,Nyb-1 do BEGIN
   Debx = i*SizeBlock
   Deby = j*SizeBlock
   Endx = (i+1)*SizeBlock-1
   Endy = (j+1)*SizeBlock-1
   if EndX GE Nx-1 then begin
        Endx = Nx-1
	Debx = EndX - SizeBlock + 1
	end
   if EndY GE Ny-1 then begin
        Endy = Ny-1
	Deby = EndY - SizeBlock + 1
	end
   TMean[ind] = (mean(Image[Debx:Endx, Deby:Endy]))
   if not keyword_set(Getnoise) then TVar[ind] = (sigma(Image[Debx:Endx, Deby:Endy]))  $
   else TVar[ind] = get_noise(Image[Debx:Endx, Deby:Endy])
   ImaSigma(i,j) = TVar[ind]
   ind = ind + 1
   END
   END
   
   ImaLarge = rebin(ImaSigma, Nxb*SizeBlock, Nyb*SizeBlock)
   ImaNoise = ImaLarge (0:nx-1,0:Ny-1)
   DataNorm = Image  /  ImaNoise 

if not keyword_set(novar) then  TVar = TVar^2
xtit='TabMean'
if not keyword_set(novar) then  ytit='TabVar' else  ytit='TabSigma'
if keyword_set(log) then begin
    TMean = alog(TMean)
    TVar = alog(TVar)
    xtit='Log TabMean'
    if not keyword_set(novar) then ytit='Log TabVar' else  ytit='Log TabSigma'
end

ind = sort(TMean)
TabVar = TVar[ind]
TabMean = TMean[ind]

FileName='fig_'+tit+'.ps'

if keyword_set(ps) then setps, filename=FileName, /portrait

if keyword_set(plot) then BEGIN
   ;window, 1
   plot, TabMean, TabVar, xtitle=xtit, ytitle=ytit, /ynozero, psym=3,title=tit, yrange=yrange, xrange=xrange
   pause, 4
   indsort = sort(TabMean)
   n10p100 = Nxb*Nyb*10/100
   Tm = (TabMean[indsort])[0:n10p100-1]
   Tv = (TabVar[indsort])[0:n10p100-1]
   ;window, 2
   ;plot, Tm, Tv, xtitle='TabMean', ytitle='TabVar', psym=3, /ynozero
   END

if not keyword_set(nofit) then begin
   res = LINFIT( TabMean, TabVar)
   a = res(1)
   a0=a
   b =  res(0)
   b0=b
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


