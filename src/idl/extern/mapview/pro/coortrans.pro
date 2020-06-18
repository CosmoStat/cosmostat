;+
; NAME:
;     COORTRANS
;
; PURPOSE:
;     Transforms between various J2000 coordinate systems.
;     No precession is performed.
;     
; CALLING SEQUENCE:
;     coortrans,coor_in,coor_out,code [, /lonlat ]
;
; INPUTS:
;     coor_in - fltarr OR dblarr - The input coordinates, specified either
;                         as (lon,lat) pairs or unit vector arrays. 
;                         Unit vectors must be of dimension (N,3).
;                         Lonlat pairs must be of dimension (N,2).
;
;     code - char string - Character string specifying the desired
;                         coordinate transformation.  Only ONE code
;                         may be requested at a time.  Valid codes are:
;                         'c2e' - in = celestial   out = ecliptic
;                         'e2c' - in = ecliptic    out = celestial
;                         'g2e' - in = galactic    out = ecliptic
;                         'e2g' - in = ecliptic    out = galactic
;                         'c2g' - in = celestial   out = galactic
;                         'g2c' - in = galactic    out = celestial
;                         'u2ll'- in = unit vector out = lon,lat
;                         'll2u'- in = lon,lat     out = unit vector
;
;
; OUTPUTS:  
;     coor_out - dblarr - The output coordinates, returned as unit
;                         vectors unless /lonlat is set by the user.
;                         Unit vectors are returned with dimension (N,3).
;                         Lonlat pairs are returned with dimension (N,2).
;
; OPTIONAL INPUT KEYWORDS:
;     /lonlat  -  Set this keyword to get output coordinates
;                         returned as longitude,latitude pairs.  
;                         (Superfluous if code = 'u2ll'). 
;
; COMMON BLOCKS:
;     None.
;
; ROUTINES CALLED:
;     get_rot_matrix()
;
; EXAMPLE:
;     Transform from celestial unit vectors (x,y,z) to Galactic (lon,lat): 
;       coortrans,[[x],[y],[z]],galcoor,'c2g',/lonlat
;     
; COMMENTS:
;     The routine does some rudimentary checking to ensure the input
;     and output formats agree with the requested transformation code.
;     It also will bounce out if the unit vectors are not normalized.
;
;     Rotation matrices are computed by the routine get_rot_matrix().
;     
; MODIFICATION HISTORY:
;     initial version, J. Weiland, 04 May 1999
;
;-
;======================================================================
;
pro coortrans,coor_in,coor_out,code,lonlat=lonlat
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
if (N_params() LT 3 ) then begin
    print,'Syntax: coortrans, coor_in, coor_out, code, [/latlon]'
    return
  endif
;
; Check to see if the input dimensions are valid: either (N,3) or (N,2)
;

dim_info = size(coor_in)
if ((dim_info[0] ne 2) OR ((dim_info[2] ne 2) and (dim_info[2] ne 3))) $
 then begin
   ;
   ; reform a single unit vector or single (lon,lat) 
   ;
   if (dim_info[0] eq 1 AND dim_info[1] eq 3) then begin
        coor_in = reform(coor_in,1,3)
        dim_info = size(coor_in)
   endif else if (dim_info[0] eq 1 AND dim_info[1] eq 2) then begin
        coor_in = reform(coor_in,1,2)
        dim_info = size(coor_in)
   endif else begin 
      print,'input coordinate dimensions must be either'
      print,' [N,2] for (longitude, latitude) pairs or [N,3] for unit vectors!'
      return
   endelse
endif
if ((dim_info[0] eq 2) AND (dim_info[2] eq 2)) then infmt = 'lonlat'
if ((dim_info[0] eq 2) AND (dim_info[2] eq 3)) then infmt = 'uvec'

;
; Ensure that input unit vectors are normalized to unit magnitude
; This could be an annoying check, since it takes up some time
; and assumes a certain accuracy for the normalization.  But leave it
; in for now.
;
if infmt eq 'uvec' then begin
  mag = sqrt(coor_in[*,0]^2 + coor_in[*,1]^2 + coor_in[*,2]^2)
  bad = where(mag gt 1.0001d0 OR mag lt 0.9998d0)
  if bad[0] ne -1 then begin
     print,'Are you sure the unit vectors are normalized?'
     return
  endif
  mag = 0           ;save space
endif
   
;
; Ensure we have a valid transformation code
;
code = strlowcase(strtrim(code,2))
if (code ne 'c2e') and (code ne 'e2c') and $
   (code ne 'c2g') and (code ne 'g2c') and $
   (code ne 'g2e') and (code ne 'e2g') and $
   (code ne 'll2u') and (code ne 'u2ll') then begin

  print,'Invalid code string: ', code
  print,'A valid code is ONE of: e2c, c2e, c2g, g2c, g2e, e2g, u2ll, ll2u '
  return
endif
;
;
; Set conversion constant between degrees and radians
;

   d2r = 3.14159265358979d0/180.d0

;
if ((infmt eq 'lonlat') or (code eq 'll2u')) then begin
   ;
   ; convert to unit vector
   ;
   lon = coor_in[*,0]
   lat = coor_in[*,1]
   uvec_in = [[cos(lon*d2r)*cos(lat*d2r)],$
              [sin(lon*d2r)*cos(lat*d2r)],$
              [sin(lat*d2r)]]

endif else begin
   uvec_in = coor_in
endelse
;
;
; Perform selected transformation
;

if (code eq 'll2u') then begin
   if infmt eq 'lonlat' then   coor_out = uvec_in
   if infmt eq 'uvec' then print,'Input coordinates are not (lon,lat)!'
   return
endif else if (code eq 'u2ll') then begin
   if infmt eq 'uvec' then begin
     ;
     ; catch case where u = [0,0,1]
     ;
     sel = where(coor_in[*,0] eq 0 AND coor_in[*,1] eq 0)
     nosel= where(coor_in[*,0] ne 0 OR coor_in[*,1] ne 0)
     coor_out = dblarr(dim_info[1],2)
     if (sel[0] ne -1) then $
       coor_out[sel,*] = [[replicate(0,n_elements(sel))],$
                          [replicate(90.,n_elements(sel))]]
     if (nosel[0] ne -1) then $
       coor_out[nosel,*] = [[atan(coor_in[*,1],coor_in[*,0])],$
                            [asin(coor_in[*,2])]]/d2r
   endif
   if infmt eq 'lonlat' then print,'Input coordinates are not unit vectors!'
   return
endif else begin
   rot_matrix = transpose(get_rot_matrix(code))
   uvec_out = rot_matrix#transpose(uvec_in)
   uvec_out = transpose(uvec_out)
endelse
;
;
; Convert to output in lonlat, if requested
;
if keyword_set(lonlat) then begin
     ;
     ; catch case where u = [0,0,1]
     ;
     sel = where(uvec_out[*,0] eq 0 AND uvec_out[*,1] eq 0)
     nosel= where(uvec_out[*,0] ne 0 AND uvec_out[*,1] ne 0)
     coor_out = dblarr(dim_info[1],2)
     if (sel[0] ne -1) then $
       coor_out[sel,*] = [[replicate(0,n_elements(sel))],$
                          [replicate(90.,n_elements(sel))]]
     if (nosel[0] ne -1) then $
       coor_out[nosel,*] = [[atan(uvec_out[*,1],uvec_out[*,0])],$
                            [asin(uvec_out[*,2])]]/d2r
   return
endif else begin
   coor_out = uvec_out
   return
endelse
;    
;
return
end






