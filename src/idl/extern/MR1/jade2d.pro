;-------------------------------------------------------------------------------
;+
; NAME:
;	jade2d
;	
; PURPOSE:
;	Apply the Independant Componant Method to a set of images using
;       the JADE method. The routine is an extensopn to 2d of jade1d routine.
;       The input set of images must be given with the following syntax:
;                       ObservSig(*, *, i) = ith image
;
; EXPLANATION:
;	compute ICA
;
; CALLING SEQUENCE:
;	jade2d, ObservSig, NbSource, DeMixingMat, Process, wave=wave, optw=optw
;
; INPUTS:
;	ObservSig : mixing of input signal (ObservSig = A # input signals)
;                    ObservSig(*, *, i) = ith vector
;       NbSource  : number of sources in input signal
;
; OPTIONAL INPUTS:
;	DispDebug : debug trace on screen
;	FileDebug : debug trace on file (Trace.dat in local directory)
;	Verbose	: verbose
;
; KEYWORD PARAMETERS:
;       wave: if set, an bi-orthogonal wavelet transform is applied
;             on each image, and the component separation is performed
;             on the wavelet transformed images.
;             The inverse wavelet transform is then applied on the
;             reconstructed process.
;       optw: string = options for the wavelet transformation
;
; OUTPUTS:
;	DeMixingMat : demixing matrix 
;       Process: Reconstructed process = DeMixingMat # ObservSig
;                Process(*,*,i) = ith Process, with i = 0.. NbSource - 1
;
; DEPENDENCIES:
;
; RESTRICTIONS:
;
; PROCEDURES USED:
;	whitening_signal, estimation_cumulant, init_diagonalisation
;	contrast_optimisation, separating_matrix
;
; EXAMPLES:
;	jade2d, ObservIma, NbSource, DeMixingMat, RecIma
;	jade2d, ObservIma, NbSource, DeMixingMat, RecIma, /wave, optw='-n5'
;
;-------------------------------------------------------------------------------   

pro jade2d, TabIma, Nprocess, Mat, TabObj, wave=wave, optw=optw, fft=fft

if N_PARAMS() LT 3 then begin 
        print, 'CALL SEQUENCE:  jade2d, ObservIma, NbSource, Ma, RecProcess, wave=wave, optw=optw'
        goto, DONE
        end
	
vs = size(TabIma)
Nx = vs[1]
Ny = vs[2]
Nz = vs[3]
Np = Nx*Ny
tab = dblarr(Nz, Np)
if keyword_set(optw) then WO = optw+' -t14 -L' $
else WO = ' -t14 -L'
MRFILE = 'xxjade.mr'
RECFILE= 'xxjade.fits'

for i=0, Nz-1 do begin
    b = TabIma(*,*,i)
    if keyword_set(wave) then begin
           ; print, WO
           mr_transform, b, w, opt=WO, MR_File_Name=MRFILE
	   b = w
	   end
    tab(i,*) = b - mean(b)
  end
; help, tab
jade, tab, Nprocess, Mat
O = Mat # tab

TabObj = dblarr(Nx, Ny, Nprocess)
if keyword_set(wave) then x = readfits(MRFILE,HD)

b = dblarr(Nx,Ny)
for i=0, Nprocess-1 do begin
   b(*) =  O(i,*)
  if keyword_set(wave) then begin
       ; help, b, HD
       writefits, MRFILE, b, HD
       COM = 'mr_recons ' + MRFILE + ' ' +  RECFILE
       spawn, com
       TabObj(*,*,i) = readfits(RECFILE)
       end $
   else TabObj(*,*,i) = b 
end

if keyword_set(wave) then begin
delete, MRFILE
delete, RECFILE
end

DONE:

end
