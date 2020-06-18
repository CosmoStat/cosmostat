; NAME:FITS2DAT
;
;
; PURPOSE: to change a .dat file (from matlab save) in a .fits file
;
;
; CALLING: 
; 					fits2dat, file_fits=file_fits, filefits, file_dat=file_dat
; INPUTS: 
;	
;	
; OUTPUT: 
;	
;
; HISTORY:
;	Written: Sandrine Pires Jan 2007.
;-
;-------------------------------------------------------------------------------

pro fits2dat, file_fits=file_fits, filefits, file_dat=file_dat

tab = readfits(file_fits)
openw, 1, file_dat
printf, 1, double(tab)
close, 1
free_lun, 1

end

