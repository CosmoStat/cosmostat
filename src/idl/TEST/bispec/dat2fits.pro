
;===================================================================================================
;
; NAME: DAT2FITS
;
;
; PURPOSE: to change a .dat file (from matlab save) in a .fits file
;
;
; CALLING: dat2fits, file_dat=file_dat filefits file_fits=file_file
;
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

pro dat2fits, file_dat=file_dat, filefits, file_fits=file_fits

tab = read_ascii(file_dat)

filefits = tab.(0)
 
writefits, file_fits, filefits

end

