;+
; NAME: 
; 				POLAR_FFT
;
;
; PURPOSE: 
; 				process the polar fft of an image
;
; CALLING: 
; 				polar_fft, f, real, ima, PolarFFT_CutoffParam=PolarFFT_CutoffParam, Name=Name
;
; INPUTS: 
;				f --- input image
;	
; OUTPUT: 
; 				real --- real part of the polar fft
; 				ima --- imaginary part of the polar fft 
;					(each line of the polar fft is a radius, 
; 				the line 0 corresponds to angle -!pi/2)
; KEYWORD:
;				PolarFFT_CutoffParam --- cut-off parameter in the approximation
; 				(the larger it is the better approximation is)
;                               Name: Name output file
;
; HISTORY:
;	Written: Sandrine Pires Sept 2009.
;-
;-------------------------------------------------------------------------------
pro polar_fft, f, real, ima, PolarFFT_CutoffParam=PolarFFT_CutoffParam, Name=Name

if keyword_set(PolarFFT_CutoffParam) then m= PolarFFT_CutoffParam else m = 8
if not keyword_set(Name) then Name='tmp'
;a = systime(1)
sz =size(f)
N = sz[1]
T=3*N; %T=5*N/2
R=3*N/2

fr=float(f) 
fi=imaginary(f)

IN_FITSFILE_RE = gettmpfilename()
IN_FITSFILE_IM = gettmpfilename()
IN_DATFILE_RE = Name+'_input_r.dat'
IN_DATFILE_IM = Name+'_input_i.dat'
OUT_DATFILE_RE = Name+'_polarfft_r.dat'
OUT_DATFILE_IM = Name+'_polarfft_i.dat'
OUT_FITSFILE_RE = gettmpfilename()
OUT_FITSFILE_IM = gettmpfilename()



writefits, IN_FITSFILE_RE, fr
fits2dat, file_fits=IN_FITSFILE_RE, filefits, file_dat=IN_DATFILE_RE
writefits, IN_FITSFILE_IM, fi
fits2dat, file_fits=IN_FITSFILE_IM, filefits, file_dat=IN_DATFILE_IM 

com = 'cea_polar_fft' +' ' +string(N)+' ' +string(T)+' ' +string(R)+' ' +string(m)+' ' +string(Name)
spawn, com

dat2fits, file_dat=OUT_DATFILE_RE, filefits, file_fits=OUT_FITSFILE_RE
dat2fits, file_dat=OUT_DATFILE_IM, filefits, file_fits=OUT_FITSFILE_IM

real = readfits(OUT_FITSFILE_RE)
ima = readfits(OUT_FITSFILE_IM)

delete, IN_FITSFILE_RE
delete, IN_FITSFILE_IM
delete, IN_DATFILE_RE
delete, IN_DATFILE_IM
delete, OUT_DATFILE_RE
delete, OUT_DATFILE_IM
delete, OUT_FITSFILE_RE
delete, OUT_FITSFILE_IM

;print, 'time =', systime(1)-a
end
