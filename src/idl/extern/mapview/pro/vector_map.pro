function vector_map, norm, Obs, vector_scale = vector_scale

;
; given an array of vector norms, 
; returns the average norm and divide norms by this average
;
;
sz = size(norm)
npix = n_elements(norm)
normout = make_array(/float, npix, value=-1.)

if defined(Obs) then begin
    if Obs[0] eq -1 then begin
        vector_scale = 0.
        return, normout
    endif
    N_obs = n_elements(Obs)
    N_no_Obs = npix - N_obs
endif else begin
    N_obs = npix
    N_no_obs = 0
endelse

if (N_NO_Obs eq 0) then begin
    if undefined(vector_scale) then begin
;     mean1 = total(norm)/N_obs
        mean1 = median(norm)
        vector_scale = mean1
    endif
    normout      = norm/vector_scale
endif else begin
    if undefined(vector_scale) then begin
;    mean1 = total(norm[Obs])/N_obs
        mean1 = median(norm[Obs])
        vector_scale = mean1
    endif
    normout[Obs] = norm[Obs]/vector_scale
endelse


return, normout
end

