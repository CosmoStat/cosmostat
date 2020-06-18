

PRO BLOCK_TRANS, Image, TabVar=TabVar, TabMean=TabMean, Plot=Plot, BlockSize=BlockSize
;+ 
; NAME: 
;       BLOCK_TRANS
;
; PURPOSE: 
;      
;       This routine separates the image in blocks of size BlockSizexBlockSize
;       and computes inside each block the mean and the variance. The plot
;       variance versus mean gives information about the kind of noise inside
;       the data. 
;
; CALLING SEQUENCE: 
;   BLOCK_TRANS, Image, TabVar=TabVar, TabMean=TabMean, 
;                Plot=Plot, BlockSize=BlockSize
;
; INPUTS: 
;   Image
;
; OPTIONAL INPUT PARAMETERS: 
;   none
;
; KEYED OUTPUT: 
;   TabVar -- IDL 1D array: array of variance
;   TabMean -- IDL 1D array: array of mean
;   Plot -- scalar: if set, the plot is displayed
;   BlockSize -- scalar: size of the blocks (default is 8)
;
; OUTPUTS: 
;   none
;
; MODIFICATION HISTORY: 
;    19-Sep-1996 JL Starck written with template_gen 
;    13-Oct-1998 bracket instead of parenthesis
;-
 
;------------------------------------------------------------
; parameters check
;------------------------------------------------------------
 
 IF N_PARAMS() LT 1 THEN BEGIN
   PRINT, 'CALLING SEQUENCE: ', $ 
    'BLOCK_TRANS, Image, TabVar=TabVar, TabMean=TabMean, Plot=Plot, BlockSize=BlockSize'
   GOTO, CLOSING
 ENDIF
 
vsize = size(Image)
if vsize[0] NE 2 then BEGIN
   print, "Error: bad first parameter ..."
   print, "       Image has to be a 2 dimensional array"
   GOTO, CLOSING
 ENDIF
Nx = vsize[1]
Ny = vsize[2]

if keyword_set(BlockSize) then SizeBlock = BlockSize else SizeBlock = 8

Nxb = Nx / SizeBlock
Nyb = Ny / SizeBlock

;------------------------------------------------------------
; function body
;------------------------------------------------------------

TMean = fltarr(Nxb*Nyb)
TVar = fltarr(Nxb*Nyb)

ind = 0
for i=0,Nxb-1 do BEGIN
for j=0,Nyb-1 do BEGIN
   Debx = i*SizeBlock
   Deby = j*SizeBlock
   Endx = (i+1)*SizeBlock-1
   Endy = (j+1)*SizeBlock-1
   TMean[ind] = (mean(Image[Debx:Endx, Deby:Endy]))
   TVar[ind] = (sigma(Image[Debx:Endx, Deby:Endy]))^2
   ind = ind + 1
   END
   END
 
ind = sort(TMean)
TabVar = TVar[ind]
TabMean = TMean[ind]

if keyword_set(plot) then BEGIN
   plot, TabMean, TabVar, xtitle='TabMean', ytitle='TabVar', /ynozero, psym=3
   pause, 4
   indsort = sort(TabMean)
   n10p100 = Nxb*Nyb*10/100
   Tm = (TabMean[indsort])[0:n10p100-1]
   Tv = (TabVar[indsort])[0:n10p100-1]
   plot, Tm, Tv, xtitle='TabMean', ytitle='TabVar', psym=3, /ynozero
   END

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


