;+
; NAME: 
;       CUR_FTRANS
;
; PURPOSE:
;       Computes the fast curvelet transform of an image. A multiresolution file
;       has a .fct extension, if the parameter file name have not one, then
;       the extension is added to the file name. Result is store in DataTransf.
;       If the keyword is not set, the created file is deleted.
;
; CALLING:
;     cur_ftrans, Imag, Trans, Opt=Opt 
;
; INPUTS:
;     Imag -- IDL 2D array: Input image be transformed 
;     
; OUTPUTS:
;     Trans -- IDL structures with the following fields:  
;         HEADTRANS -- IDL STRING Array: contains the Fits Header of the decomposition
;
; KEYWORDS:
;
;      Opt -- string: string which contains the differents options. Options are:
;
;        where options = 
;
;   
;           [-n number_of_scales]
;                number of scales used in the wavelet transform
;                default is 4
;
;           [-r]
;                Inverse Fast Curvelet transform.
;
;           [-v]
;                Verbose. Default is no.
;
; EXTERNAL CALLS:
;       cur_ftrans (C++ program)
;
; EXAMPLE:
;
; HISTORY:
;       Written : Jean-Luc Starck, October 2005.
;-
;-----------------------------------------------------------------

; b = band in  [1..Cur.NbrBand]
; s = scale in [1..Cur.NbrScale]
; d = direction [1.. Cur.TabNbrAnglePerScale[s-1]

function fctget, Cur, b=b, s=s, d=d, re=re, im=im, pow=pow
  if keyword_set(s) then begin
     if s LT 1 or s GT Cur.NbrScale then begin
         print, 'Error: Scale must be between 1 and ', Cur.NbrScale
        goto, DONE
        end
  end
  if keyword_set(b) then begin
     if b LT 1 or b GT Cur.NbrBand then begin
         print, 'Error: Band must be between 1 and ', Cur.NbrBand
        goto, DONE
        end
  end
  if keyword_set(d) then begin
     if not keyword_set(s) then begin
        print, 'Error: Scale keyword must be set'
        goto, DONE
        end
    
     if d LT 1 or d GT Cur.TabNbrAnglePerScale[s-1] then begin
         print, 'Error: Direction number must be between 1 and ', Cur.TabNbrAnglePerScale[s-1]
        goto, DONE
        end
  end
  	
  FirstPos=11
  IndBand=FirstPos
  if keyword_set(b) then begin
    IndBand=FirstPos+b-1
    BandNumber = b-1 
  end else if  keyword_set(s) then begin
     Nb=0
     for i=0,s-2 do  Nb=Nb+Cur.TabNbrAnglePerScale[i]
     IndBand=IndBand+Nb
     if keyword_set(d) then IndBand=IndBand+d-1
     BandNumber = IndBand - FirstPos
  end
  ; print, s, d, IndBand
  band = cur.(IndBand)
  ; help, band
  b=BandNumber
  ; print, cur.tabbandnx(b), cur.tabbandny(b), cur.tabbandnx(b)*cur.tabbandny(b)
  ; band = reform(Band, cur.tabbandnx(b), cur.tabbandny(b))
  if Cur.Real NE 1 then begin
     if keyword_set(re) then band = float(band) 
     if keyword_set(im) then band = imaginary(band)
     if keyword_set(pow) then band = (float(band))^2 + (imaginary(band))^2
  end
  
  return, band
DONE:

end

;====================================


pro fctput, Cur, Band, b=b, s=s, d=d, re=re, im=im 
  if keyword_set(s) then begin
     if s LT 1 or s GT Cur.NbrScale then begin
         print, 'Error: Scale must be between 1 and ', Cur.NbrScale
        goto, DONE
        end
  end
  if keyword_set(b) then begin
     if b LT 1 or b GT Cur.NbrBand then begin
         print, 'Error: Band must be between 1 and ', Cur.NbrBand
        goto, DONE
        end
  end
  if keyword_set(d) then begin
     if not keyword_set(s) then begin
        print, 'Error: Scale keyword must be set'
        goto, DONE
        end
    
     if d LT 1 or d GT Cur.TabNbrAnglePerScale[s-1] then begin
         print, 'Error: Direction number must be between 1 and ', Cur.TabNbrAnglePerScale[s-1]
        goto, DONE
        end
  end
  	
  FirstPos=11
  IndBand=FirstPos
  if keyword_set(b) then IndBand=FirstPos+b-1 $
  else if  keyword_set(s) then begin
     Nb=0
     for i=0,s-2 do  Nb=Nb+Cur.TabNbrAnglePerScale[i]
     IndBand=IndBand+Nb
     if keyword_set(d) then IndBand=IndBand+d-1
  end
  ; print, s, d, IndBand
  
  cur.(IndBand) = band ; reform(band, cur.tabbandnx(indband)*cur.tabbandny(indband))
  
  if cur.real NE 1 then BEGIN
    if keyword_set(im) then begin
       bandRe = float(cur.(IndBand))
       cur.(IndBand) = complex(bandRe, Band)
       end else if keyword_set(re) then begin
          bandIm = imaginary(band)
          cur.(IndBand) = complex(band, BandIm)
       end else cur.(IndBand) = band
  END
DONE:

end

;====================================

pro test_curstat, Cur
for s=1, Cur.NbrScale do begin
for d=1, Cur.TabNbrAnglePerScale[s-1] do begin
   print, 'Scale', s, ' Dir ', d, sigma(fctget(Cur, s=s, d=d, /re)), sigma(fctget(Cur, s=s, d=d, /im))
end
end
end

;====================================


pro cur_frec, CurTrans, Imag

Header=CurTrans.HeadTrans
Cur=fltarr(CurTrans.NxCur)
Cur[0] = CurTrans.NbrScale
Cur[1:CurTrans.NbrScale] = CurTrans.TabNbrAnglePerScale  
ind = 1L + CurTrans.NbrScale
indband = 0
Ns = CurTrans.NbrScale
for s=0,CurTrans.NbrScale-1 do begin
for b=0,CurTrans.TabNbrAnglePerScale[s]-1 do begin
  Ny = CurTrans.TabBandNy[indband]
  Nx = CurTrans.TabBandNx[indband]   
  Cur[ind] = Ny
  Cur[ind+1] = Nx
  if CurTrans.real NE 1 or s EQ Ns-1 then a = fctget(CurTrans, s=s+1, d=b+1, /re) $
  else  a = fctget(CurTrans, s=s+1, d=b+1)
  ;print, s, b, Nx, Ny
  ;help, a
  Cur[ind+2:ind+1+Nx*Ny] = reform(a, Nx*Ny)
  ind = ind + 2 + Nx*Ny
  if CurTrans.real NE 1 then BEGIN
    Cur[ind:ind+Nx*Ny-1] = reform(fctget(CurTrans, s=s+1, d=b+1, /im), Nx*Ny)
    ind = ind + Nx*Ny
  END
  indband = indband + 1
end
end

; filename = 'xx_mr.fct'
filename = 'xx.mr'
NameImag='xx_imag.fits'
writefits, filename, Cur, Header
; com = 'cur_ftrans -r '  + filename + ' ' +  NameImag
com = 'mr_recons '  + filename + ' ' +  NameImag
; print, com
spawn, com
Imag = readfits(NameImag)

; delete, NameImag
; delete, filename
DONE:
end


;====================================

pro cur_ftrans, Ima, CurTrans, OPT=OPT, real=real

if N_PARAMS() LT 2 then begin 
        spawn, 'cur_trans'
        print, 'CALL SEQUENCE: cur_trans, Signal, Struct_Out, OPT=Opt'
        goto, DONE
        end

Nx = (size(Ima))[1]
Ny = (size(Ima))[2]

NameIma = 'xx_signal.fits'
NameResult = 'xx_result'
NameResultFits = 'xx_result.fct'

writefits,  NameIma,  Ima
if not keyword_set(OPT) then OPT = ' '
; if   keyword_set(real) then OPT = OPT  + ' -R '

;com = 'cur_ftrans'+' '+ OPT + ' '+ NameIma  + ' ' +  NameResult
NameResultFits = 'xx_result.mr'
com = 'mr_transform -t28 -v '+' '+ OPT + ' '+ NameIma  + ' ' +  NameResult

; print, com
spawn, com

Cur = readfits(NameResultFits, HeadTrans, /silent);         
Nl = FXPAR(HeadTrans, "NL")
NC = FXPAR(HeadTrans, "NC")
NBRDIR = FXPAR( HeadTrans, "NBRDIR")
NbrScale = FXPAR( HeadTrans, "NbrScale")
Real = FXPAR( HeadTrans, "Real")
NxCur = (size(Cur))[1]
Ns = long(Cur[0])
TabNbrAnglePerScale = long(Cur[1:Ns])
NbrBand = total( TabNbrAnglePerScale )
TabDepRe = lonarr(NbrBand)
TabEndRe = lonarr(NbrBand)
TabDepIm = lonarr(NbrBand)
TabEndIm = lonarr(NbrBand)
TabBandNx = lonarr(NbrBand)
TabBandNy = lonarr(NbrBand)
BandName = strarr(NbrBand)
ind = 1L + Ns
indband = 0L
for s=0L,Ns-1 do begin
for b=0L,TabNbrAnglePerScale[s]-1 do begin
  Ny = long(Cur[ind])
  Nx = long(Cur[ind+1])
  Dep = ind + 2L
  TabBandNx[indband] = Nx
  TabBandNy[indband] = Ny
  TabDepRe[indband] = Dep
  TabEndRe[indband] = Dep+Nx*Ny-1L
  C1 =Cur[ TabDepRe[indband]:TabEndRe[indband]]
  ; print, Nx, Ny, nx*ny, Dep, Dep+Nx*Ny-1L - Dep + 1
  BandRe =  reform(C1, TabBandNx[indband], TabBandNy[indband])
  ind = ind + 2L + Nx*Ny
  if Real NE 1 then begin
    TabDepIm[indband] = ind
    TabEndIm[indband] = ind+Nx*Ny-1
    ind = ind + Nx*Ny 
    C1 = Cur[ TabDepIm[indband]:TabEndIm[indband]]
    BandIm =  reform(C1, TabBandNx[indband], TabBandNy[indband])
  end
  my_command = 'scale_'+strcompress(string(s+1), /remove_all) + '_dir_' + strcompress(string(b+1), /remove_all)
  BandName[indband] = my_command
  ; my_command = BandName[indband] + ' = complex(Cur[ TabDepRe[indband]:TabEndRe[indband]], Cur[ TabDepIm[indband]:TabEndIm[indband]])'
  if Real NE 1 then my_command = BandName[indband] + ' = complex(BandRe, BandIm)' $
  else my_command = BandName[indband] + ' =  BandRe'
  ; print,  TabDepRe[indband], TabEndRe[indband], TabDepIm[indband], TabEndIm[indband]
  ; print, my_command
  ACK = EXECUTE( my_command)
  indband = indband + 1L 
 end
end
; help, scale_3_dir_6

;indband = 0
;ind = 1L + Ns
;for s=0,Ns-1 do begin
;for b=0,TabNbrAnglePerScale[s]-1 do begin
;  Ny = Cur[ind]
;  Nx = Cur[ind+1] 
;  BandR = fltarr(Nx,Ny)
;  BandR[*] = Cur[ind+2:ind+1+Nx*Ny]
;  ind = ind + 2 + Nx*Ny
;  BandI = fltarr(Nx,Ny)
;  BandI[*] = Cur[ind:ind+Nx*Ny-1]
;  print, sigma(BandR), sigma(BandI), sigma(  Cur[ TabDepRe[indband]:TabEndRe[indband]] ), sigma(  Cur[ TabDepIm[indband]:TabEndIm[indband]] )
;  ind = ind + Nx*Ny
;  indband = indband + 1
;end
;end
 
  my_command = 'CurTrans = { NbrScale : NbrScale,'
  my_command = my_command +'NBRDIR : NBRDIR, '
  my_command = my_command +'TabNbrAnglePerScale : TabNbrAnglePerScale, ' 
  my_command = my_command +'NbrBand : NbrBand, '
  my_command = my_command +'Real : Real, '
  my_command = my_command +'TabDepRe: TabDepRe, '
  my_command = my_command +'TabEndRe: TabEndRe, '
  my_command = my_command +'TabDepIm: TabDepIm, '
  my_command = my_command +'TabEndIm: TabEndIm, '
  my_command = my_command +'TabBandNx: TabBandNx, '
  my_command = my_command +'TabBandNy: TabBandNy, '
  indband = 0
  for s=0,Ns-1 do begin
  for b=0,TabNbrAnglePerScale[s]-1 do begin
    my_command = my_command + BandName[indband] +':'+ BandName[indband] +','
    indband = indband + 1
  end
  end
  my_command = my_command + ' HeadTrans:HeadTrans,'
  my_command = my_command + ' NxCur :NxCur'
  my_command = strcompress( my_command, /remove_all)
  my_command =my_command+'}'
  ; print, my_command
  ACK = EXECUTE( my_command)
  
  
;CurTrans = {Nxcur : Nxcur, Nycur:Nycur, Nzcur:Nzcur,Coef : Cur, Bsize : Bsize, $
;            Overlap: Overlap, HeadTrans:HeadTrans}
DONE:
 
end


;====================================


pro cur_fstat, CurTrans, Imag
 
ind = 1L + CurTrans.NbrScale
indband = 0
print, CurTrans.NbrBand
TabSigRe = fltarr(long(CurTrans.NbrBand))
TabSigIm = fltarr(long(CurTrans.NbrBand))
help, TabSigRe
for s=0,CurTrans.NbrScale-1 do begin
for b=0,CurTrans.TabNbrAnglePerScale[s]-1 do begin
  Ny = CurTrans.TabBandNy[indband]
  Nx = CurTrans.TabBandNx[indband]   
  bandr = fctget(CurTrans, s=s+1, d=b+1, /re)
  TabSigRe[indband]  = sigma(bandr)
  M1 = mean(bandr)
  if CurTrans.Real NE 1 then BEGIN
     Bandi = fctget(CurTrans, s=s+1, d=b+1, /im)
     TabSigIm[indband]  = sigma(bandi)
     M2 = mean(bandi)
     print, "BAND " , indband,  TabSigRe[indband], TabSigIm[indband], " ", M1, " ", M2
  END else  print, "BAND " , indband,  TabSigRe[indband], " ", M1 
  indband = indband + 1
end
end
m = max(TabSigRe, pos)
print, "BAND RE MAX = ", pos, m 
DONE:
end

;====================================

pro cur_scale_fstat, CurTrans, Scale, TabSigRe, TabSigIm, AllBandRe=AllBandRe, AllBandIm=AllBandIm, verb=verb
ind = 1L + CurTrans.NbrScale
indband = 0
if keyword_set(verb) then print, 'Nbr Band = ', CurTrans.NbrBand
TabSigRe = fltarr(long(CurTrans.TabNbrAnglePerScale[Scale]))
TabSigIm = fltarr(long(CurTrans.TabNbrAnglePerScale[Scale]))
s=Scale
for b=0,CurTrans.TabNbrAnglePerScale[s]-1 do begin
  Ny = CurTrans.TabBandNy[indband]
  Nx = CurTrans.TabBandNx[indband]   
  bandr = fctget(CurTrans, s=s+1, d=b+1, /re)
  TabSigRe[indband]  = sigma(bandr)
  if CurTrans.Real NE 1 then BEGIN
     Bandi = fctget(CurTrans, s=s+1, d=b+1, /im)
     TabSigIm[indband]  = sigma(bandi)
  END
  if keyword_set(verb) then print, "   BAND RE" , indband,  'Sigma = ', TabSigRe[indband], ' Skewness = ', skewness(bandr), ' Kurtosis = ', kurtosis(bandr)  
  if (keyword_set(verb) and CurTrans.Real NE 1)   then print, "   BAND IM" , indband,  'Sigma = ', TabSigIm[indband], ' Skewness = ', skewness(bandi), ' Kurtosis = ', kurtosis(bandi) 
  indband = indband + 1
  if b EQ 0 then AllBandRe = reform(bandr, N_ELEMENTS(bandr)) else AllBandRe = [AllBandRe, reform(bandr, N_ELEMENTS(bandr))]
  if (b EQ 0 and CurTrans.Real NE 1) then AllBandIm = reform(bandi, N_ELEMENTS(bandi)) else AllBandIm = [AllBandIM, reform(bandi, N_ELEMENTS(bandi))]
end

DONE:

end

;====================================

pro cur_allscale_fstat, CurTrans

for j=0, CurTrans.NbrScale-1 do begin
cur_scale_fstat, CurTrans, j, TabSigRe, TabSigIm, AllBandRe=bandr, AllBandIm=bandi
print, "   Scale RE" , j+1,  '     Sigma = ', sigma(bandr), ' Skewness = ', skewness(bandr), ' Kurtosis = ', kurtosis(bandr)  
if CurTrans.Real NE 1 then print, "   Scale IM" , j+1,  '     Sigma = ', sigma(bandi), ' Skewness = ', skewness(bandi), ' Kurtosis = ', kurtosis(bandi) 
end

end

;====================================

pro testrec, n, opt=opt
; n = rim('lena_li.fits')

; n = rim('ngc2997.fits')

; n = rim('n.fits')
cur_ftrans, n, c, /re, opt=opt  
cur_frec, c, r
info, n-r
;info,n -t
 load, r
;cb = readfits('tt.fct')
;cc = readfits('xx_mr.fct')
end

