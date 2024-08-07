;+
; NAME: 
;       WT2D1D_Trans
;
; PURPOSE:
;       Computes the wavelet transform of a cube, by first applying a 2D WT on each
;       frame of the cube, and then a 1D WT along z. A multiresolution file
;       has a .mr extension, if the parameter file name have not one, then
;       the extension is added to the file name. Result is store in DataTransf.
;       If the keyword is not set, the created file is deleted.
;
; CALLING:
;     WT2D1D_Trans, Cube, Trans, Opt=Opt 
;
; INPUTS:
;     Cube -- IDL 3D array: Input cube be transformed 
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
;           [-n number_of_scales_2D]
;                number of scales used in the wavelet transform
;                default is 5
;           [-N number_of_scales_1D]
;                number of scales used in the wavelet transform
;                default is 4
;
;           [-r]
;                Inverse  transform.
;
;           [-v]
;                Verbose. Default is no.
;
; EXTERNAL CALLS:
;       mr2d1d_trans (C++ program)
;
; EXAMPLE:
;
; HISTORY:
;       Written : Jean-Luc Starck, January 2007.
;-
;-----------------------------------------------------------------


;====================================


function wt2d1dget, Trans, b=b, s2d=s2d, s1d=s1d
  if keyword_set(b) then begin
     if b LT 1 or b GT Trans.NbrBand then begin
         print, 'Error: Band must be between 1 and ', Trans.NbrBand
        goto, DONE
        end
  end
  if keyword_set(s2d) then begin
     if s2d LT 1 or s2d GT Trans.NbrBand2D then begin
         print, 'Error: Scale must be between 1 and ', Trans.NbrBand2D
        goto, DONE
        end
  end
  if keyword_set(s1d) then begin
     if s1d LT 1 or s1d GT Trans.NbrBand1D then begin
         print, 'Error: Band must be between 1 and ', Trans.NbrBand1D
        goto, DONE
        end
  end

  FirstPos=10
  IndBand=FirstPos
  if keyword_set(b) then begin
    IndBand=FirstPos+b-1
    BandNumber = b-1 
  end else if  keyword_set(s2d) then begin
     Nb=0
     for i=0,s2d-2 do  Nb=Nb+ Trans.NbrBand1D
     IndBand=IndBand+Nb
     if keyword_set(s1d) then IndBand=IndBand+s1d-1
     BandNumber = IndBand - FirstPos
  end
     ; print, "Ind Band GET = ", IndBand

  ; print, s, d, IndBand
  band = Trans.(IndBand)
  ; help, band
  b=BandNumber
  ; print, cur.tabbandnx(b), cur.tabbandny(b), cur.tabbandnx(b)*cur.tabbandny(b)
  ; band = reform(Band, cur.tabbandnx(b), cur.tabbandny(b))
 
  
  return, band
DONE:

end

;====================================

pro wt2d1put, Trans, Band, b=b, s2d=s2d, s1d=s1d
   if keyword_set(b) then begin
     if b LT 1 or b GT Trans.NbrBand then begin
         print, 'Error: Band must be between 1 and ', Trans.NbrBand
        goto, DONE
        end
  end
  if keyword_set(s2d) then begin
     if s2d LT 1 or s2d GT Trans.NbrBand2D then begin
         print, 'Error: Scale must be between 1 and ', Trans.NbrBand2D
        goto, DONE
        end
  end
  if keyword_set(s1d) then begin
     if s1d LT 1 or s1d GT Trans.NbrBand1D then begin
         print, 'Error: Band must be between 1 and ', Trans.NbrBand1D
        goto, DONE
        end
  end

  FirstPos=10
  IndBand=FirstPos
  if keyword_set(b) then IndBand=FirstPos+b-1 $
  else if  keyword_set(s2d) then begin
     Nb=0
     for i=0,s2d-2 do  Nb=Nb+ Trans.NbrBand1D
     IndBand=IndBand+Nb
     if keyword_set(s1d) then IndBand=IndBand+s1d-1
  end
   ; print, "Ind Band PUT= ", IndBand
  Trans.(IndBand) =  band 
DONE:

end


function strc, x
return, STRCOMPRESS(string(X), /REMOVE_ALL)
end


;====================================

pro wt2d1d_stat, Trans
for s=1, Trans.NbrBand2D do begin
for s1=1, Trans.NbrBand1D do begin
   Band = wt2d1dget(Trans, s2d=s, s1d=s1)
   vs = size(Band)
   
   sz = '(' + strc(vs[1]) + ',' + strc(vs[2]) + ',' + strc(vs[3])+ ')'
   print, 'Band2D', s, ' Band1D ', s1, ' Size = ', sz, '     Sigma = ', sigma(band), ' Skewness = ', skewness(band), ' Kurtosis = ', kurtosis(band)
end
end
end

;====================================

pro wt2d1d_rec, Trans, Imag, WT

Header=Trans.HeadTrans
NbrElem =  FXPAR(Header, "NAXIS1")

WT=fltarr(NbrElem)
WT[0] = Trans.NbrBand2D
WT[1] = Trans.NbrBand1D 
ind = 2L  
indband = 0
for s=0,Trans.NbrBand2D-1 do begin
for s1=0,Trans.NbrBand1D-1 do begin
  Ny = Trans.TabBandNy[indband]
  Nx = Trans.TabBandNx[indband]  
  Nz = Trans.TabBandNz[indband]   
  WT[ind] = Nx
  WT[ind+1] = Ny
  WT[ind+2] = Nz
  a = wt2d1dget(Trans, s2d=s+1, s1d=s1+1) 
  WT[ind+3:ind+2+Nx*Ny*Nz] = reform(a, Nx*Ny*Nz)
  ind = ind + 3 + Nx*Ny*Nz
  indband = indband + 1
end
end

filename = 'xx_mr.mr'
NameImag='xx_imag.fits'
writefits, filename, WT, Header
com = 'mr2d1d_trans -r '  + filename + ' ' +  NameImag
; print, com
spawn, com
Imag = readfits(NameImag)

delete, NameImag
delete, filename

DONE:
end


;====================================

pro wt2d1d_trans, Cube, Trans, OPT=OPT 

if N_PARAMS() LT 2 then begin 
        spawn, 'wt2d1d_trans'
        print, 'CALL SEQUENCE: wt2d1d_trans, Signal, Struct_Out, OPT=Opt'
        goto, DONE
        end

Nx = (size(Cube))[1]
Ny = (size(Cybe))[2]
Nz = (size(Cube))[3]

NameIma = 'xx_signal.fits'
NameResult = 'xx_result'
NameResultFits = 'xx_result.mr'

writefits,  NameIma,  Cube
if not keyword_set(OPT) then OPT = ' '
 
com = 'mr2d1d_trans'+' '+ OPT + ' '+ NameIma  + ' ' +  NameResult
; print, com
spawn, com

WT2D1D = readfits(NameResultFits, HeadTrans, /silent); 
delete, NameIma
delete, NameResultFits
        
Nx = FXPAR(HeadTrans, "NX")
NY = FXPAR(HeadTrans, "NY")
NZ = FXPAR(HeadTrans, "NZ")
NbrScale2D = FXPAR( HeadTrans, "NSCALE2D")
NbrScale1D = FXPAR( HeadTrans, "NSCALE1D")
Trans =  FXPAR( HeadTrans, "TYPE_TRA")
NP = (size(WT2D1D))[1]
NbrBand2D = fix(WT2D1D[0])
NbrBand1D = fix(WT2D1D[1])
NbrBand = NbrBand1D * NbrBand2D
TabDep = lonarr(NbrBand)
TabEnd = lonarr(NbrBand)
TabBandNx = lonarr(NbrBand)
TabBandNy = lonarr(NbrBand)
TabBandNz = lonarr(NbrBand)
BandName = strarr(NbrBand)
indband = 0L
Ind=2L
for s=0,NbrBand2D-1 do begin
for s1=0,NbrBand1D-1 do begin
  Nxb = long( WT2D1D[ind] )
  Nyb = long( WT2D1D[ind+1L])
  Nzb = long(WT2D1D[ind+2L]) 
  TabBandNx[indband] = Nxb
  TabBandNy[indband] = Nyb
  TabBandNz[indband] = Nzb
  SizeBand = Nxb * Nyb * Nzb
  
  TabDep[indband] = ind+3L
  TabEnd[indband] = TabDep[indband] + SizeBand -1L
  
 ; print, "ind = ", Ind, Nxb, Nyb, Nzb ;  WT2D1D[ TabEnd[indband]+1 ], WT2D1D[ TabEnd[indband]+2 ], WT2D1D[ TabEnd[indband]+3 ]
 ; print, s+1, s1+1, TabBandNx[indband], TabBandNy[indband], TabBandNz[indband]
  Coef = WT2D1D[ TabDep[indband]:TabEnd[indband]] 
;  help, Coef
;  print , SizeBand
  
  Band =  reform(Coef,  TabBandNx[indband], TabBandNy[indband], TabBandNz[indband])
  Ind = TabEnd[indband] + 1L
 
  my_command = 's2d_'+strcompress(string(s+1), /remove_all) + '_s1d_' + strcompress(string(s1+1), /remove_all)
  BandName[indband] = my_command
  ; my_command = BandName[indband] + ' = complex(WT2D1D[ TabDepRe[indband]:TabEndRe[indband]], WT2D1D[ TabDepIm[indband]:TabEndIm[indband]])'
  my_command = BandName[indband] + ' =  Band'
  ; print,  TabDepRe[indband], TabEndRe[indband], TabDepIm[indband], TabEndIm[indband]
  ; print, my_command
  ACK = EXECUTE( my_command)
  indband = indband + 1L
 end
end
 
  my_command = 'Trans = { NbrScale2D : NbrScale2D,'
  my_command = my_command +'NbrScale1D : NbrScale1D, '
  my_command = my_command +'NbrBand2D : NbrBand2D, '
  my_command = my_command +'NbrBand1D : NbrBand1D, '
  my_command = my_command +'NbrBand : NbrBand, '
  my_command = my_command +'TabDep: TabDep, '
  my_command = my_command +'TabEnd: TabEnd, '
  my_command = my_command +'TabBandNx: TabBandNx, '
  my_command = my_command +'TabBandNy: TabBandNy, '
  my_command = my_command +'TabBandNz: TabBandNz, '
  indband = 0
  for s=0,NbrBand2D-1 do begin
  for s1=0,NbrBand1D-1 do begin
    my_command = my_command + BandName[indband] +':'+ BandName[indband] +','
    indband = indband + 1
  end
  end 
  my_command = my_command + ' HeadTrans:HeadTrans,'
  my_command = my_command + ' Nx :Nx,'
  my_command = my_command + ' Ny :Ny,'
  my_command = my_command + ' Nz :Nz'
  my_command = strcompress( my_command, /remove_all)
  my_command =my_command+'}'
  ; print, my_command
  ACK = EXECUTE( my_command)
   
DONE:
 
end
 
 
;====================================

;====================================

 pro wt2d1d_test
 Cube = randomn(seed, 64,64,128)*100.
 wt2d1d_trans, Cube, Trans, OPT='-t2 -v -n3 -N4 -M'
 ; wt2d1d_stat, Trans
 
 for s=1, Trans.NbrBand2D do begin
 for s1=1, Trans.NbrBand1D do begin
    Band = wt2d1dget(Trans, s2d=s, s1d=s1)   
    wt2d1put, Trans, Band, s2d=s, s1d=s1
 end
 end
 wt2d1d_rec, Trans, rec
 info, cube-rec
 end
 
;====================================

pro wt2d1d_filter, Cube, Rec, OPT=OPT, Nsigma=Nsigma, Trans=Trans

if not keyword_set(OPT) then opt=' -t2 -A '
if not keyword_set(Nsigma) then Nsigma = 3.

; Cube = randomn(seed, 64,64,256)*100.
wt2d1d_trans, Cube, Trans, OPT=OPT
for s=1, Trans.NbrBand2D do begin
for s1=1, Trans.NbrBand1D do begin
   Band = wt2d1dget(Trans, s2d=s, s1d=s1)   
   ; print, 'Band2D', s, ' Band1D ', s1, '     Sigma = ', sigma(band), ' Skewness = ', skewness(band), ' Kurtosis = ', kurtosis(band)

   Nel = N_ELEMENTS(Band)
   SigmaNoise = mad(Band)
   ThresholdLevel = SigmaNoise * Nsigma
   index = where ( ABS(Band) LT ThresholdLevel, count )
   if count GT 0 then Band[index] = 0
   print, 'Band2D', s, ' Band1D ', s1, ": Thres = ", ThresholdLevel, ", Percentage of Thres = ", float(count) / float(Nel) * 100.
   wt2d1put, Trans, Band, s2d=s, s1d=s1
end
end

wt2d1d_rec, Trans, rec
end

;====================================

pro test_2d1dfilter

cube = readfits('test.fits')


wt2d1d_filter, Cube, Rec3, Nsigma=3, opt='-t24 '
writefits, 'rec3.fits', Rec3

disp
load, total(Cube,3)
Disp, win=1
load, total(Rec3,3)
disp, win=2
load, total(Cube-Rec3,3)


wt2d1d_filter, Cube, Rec5, Nsigma=5, opt='-t24 ' 
writefits, 'rec5.fits', Rec5


Disp, win=3
load, total(Rec5,3)
disp, win=4
load, total(Cube-Rec5,3)




end