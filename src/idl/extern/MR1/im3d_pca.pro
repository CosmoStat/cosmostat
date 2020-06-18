;-------------------------------------------------------------------------------
;+
; NAME:
;        IM3D_PCA
;
; PURPOSE:
;	Computes the Principal Componant Analysis method to a set of images.
;       The input set of images must be given with the following syntax:
;                       ObservSig(*, *, i) = ith image
;
; CALLING:
;
;      IM3D_PCA, ObservSig, PCAOut, PCA_File_Name=PCA_File_Name, RecPCA=RecPCA, NFirstPCA=NFirstPCA, wave=wave, optw=optw
;
; INPUTS:
;     ObservSig : 3D fltarr = input data. ObservSig(*, *, i) = ith vector;      
;
; OUTPUTS:
;     PCAOut: 3D fltarr = output data.  PCAOut(*, *, i) = ith eigen vector; 
;
; KEYWORDS:
;      PCA_File_Name: string containing the PCA file name
;      wave: if set, a bi-orthogonal wavelet transform is applied
;             on each image, and the PCA is performed
;             on the wavelet transformed images.
;             The inverse wavelet transform is then applied on the eigen vectors.
;      optw: string = options for the wavelet transformation
;
;      RecPCA: output 3D fltarr: reconstructed data from the NFirstPCA eigen vectors
;
;      NFirstPCA: integer: number of eigen vectors used for the reconstruction
;
; EXTERNAL CALLS:
;       mr_transform, mr_recons, im3d_pca (C++ program)
;
; EXAMPLE:
;
;       Compute the PCA of the cube. The result is stored in 
;       the file result.fits
;               IM3D_PCA, Data, PCA, PCA_File_Name='result.fits'
;
;       Compute the PCA of the cube and reconstruct a cube from
;       the first three eigen vectors.
;              IM3D_PCA,, Data, PCA,  RecPCA=RecPCA, OPT='-F 3'
;
; HISTORY:
;	Written: Jean-Luc Starck 2004.
;	July, 2004 File creation
;-

pro  IM3D_PCA, ObservSig, PCAOut, PCA_File_Name=PCA_File_Name, RecPCA=RecPCA, NFirstPCA=NFirstPCA, wave=wave, optw=optw

ImagIn = float(ObservSig)
vs = size(ImagIn)

if N_PARAMS() NE 2 then begin 
        print, 'CALLING SEQUENCE: IM3D_PCA, ObservSig, PCAOut, PCA_File_Name=PCA_File_Name, RecPCA=RecPCA, NFirstPCA=NFirstPCA, wave=wave, optw=optw'
        goto, DONE
        end

if vs[0] NE 3 then begin
        print, 'ERROR: bad first parameter ...'
        print, 'CALLING SEQUENCE: IM3D_PCA, ObservSig, PCAOut, PCA_File_Name=PCA_File_Name, RecPCA=RecPCA, NFirstPCA=NFirstPCA, wave=wave, optw=optw'
        goto, DONE   
        end
	
Nx = vs[1]
Ny = vs[2]
Nz = vs[3]
Np = Nx*Ny
tab = dblarr(Nz, Np)
if keyword_set(optw) then WO = optw+' -t14 -L' $
else WO = ' -t14 -L'
MRFILE = 'xx_pca.mr'
RECFILE= 'xx_pca.fits'

if keyword_set(wave) then begin
    for i=0, Nz-1 do begin
       b = ImagIn(*,*,i)
       mr_transform, b, w, opt=WO, MR_File_Name=MRFILE
       ImagIn(*,*,i) = w
   end
   x = readfits(MRFILE,HD)
end
  
if not keyword_set(Optw) then Optw = ' '
if keyword_set(PCA_File_Name) then NF=1 else NF = 0
if keyword_set(NFirstPCA) then NF=NFirstPCA 

if keyword_set(PCA_File_Name) then filename = PCA_File_Name else filename = 'xx_temp.fits'
PcaRec_Filename='xx_rec.fits'

if NF gt 0 then  opt='-r xx_rec.fits' +  ' -F ' + STRCOMPRESS(STRING(NF), /REMOVE_ALL) $
else opt=' '

p = strpos(filename, '.fits')
if p LT 0 then filename = filename + '.fits'

NameImag='xx_imag.fits'
writefits, NameImag, ImagIn

com = 'im3d_pca ' + ' ' + Opt + ' ' + NameImag + ' ' +  filename
spawn, com
PCAOut = readfits(filename)
if keyword_set(NFirstPCA) then RecPCA= readfits('xx_rec.fits')

if keyword_set(wave) then begin
    for i=0, Nz-1 do begin
       b = PCAOut (*,*,i)
       writefits, MRFILE, b, HD
       COM = 'mr_recons ' + MRFILE + ' ' +  RECFILE
       spawn, com
       PCAOut (*,*,i) = readfits(RECFILE)
    end
    if keyword_set(NFirstPCA) then begin
    for i=0, Nz-1 do begin
       b = RecPCA (*,*,i)
       writefits, MRFILE, b, HD
       COM = 'mr_recons ' + MRFILE + ' ' +  RECFILE
       spawn, com
       RecPCA (*,*,i) = readfits(RECFILE)
    end
    end
end
 
delete, NameImag
if filename EQ 'xx_temp.fits' then delete, filename
if keyword_set(NFirstPCA) then delete, 'xx_rec.fits'
DONE:

END
