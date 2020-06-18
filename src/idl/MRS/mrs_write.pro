;+
; NAME:
;        mrs_write
;
; PURPOSE:
;	
;   write, spherical map, either in healpix or glesp format
;
;
; CALLING:
;
;      mrs_write, file, data, ring=ring
;       
; INPUTS:
;    file : file to be writen
;   data : healpix tab or glesp struct to be writen for Healpix map the map is in NESTED format
;    
; OUTPUTS:
;     
;    none
; KEYWORDS:
;    ring: convert Healpix map data to RING format for writing
;
; EXTERNAL CALLS:
;       write_glesp (glesp)
;   	write_fits_map (healpix)
;
; EXAMPLE:
;
;   mrs_write,"my_file_healpix.fits",healpixdata
;   mrs_write,"my_file_glesp.fits",glespdata
;         
; HISTORY:
;	Written: Pierrick Abrial & Jean-Luc Starck, 2006
;	February 2005
;
;
;		mrsp_write, file, data, ring=ring	Write a HealPIX (no glesp) POLARIZED map in NESTED format with the same options and calling as mrs_write
;
;
;-------------------------------------------------------------------------------------------

pro mrsp_write, file, data, ring=ring

     if not keyword_set(ring) then  write_tqu, file,data,/nested  $       
    else begin 
      map_ring = reorder(data,out='ring',in='nest')
      write_tqu, file, map_ring,/ring
    end
    
 
END


;========================================================================

pro mrs_write, file, data, ring=ring

  if type_code(data) EQ 8 then begin ; type = glesp
    write_glesp,file,data
    endif else begin 
    
    if not keyword_set(ring) then  write_fits_map,file,data,/nested  $       
    else begin 
      map_ring = reorder(data,out='ring',in='nest')
      write_fits_map,file,map_ring,/ring
    end
 endelse 
  
END

