pro df_mask_lines, inlam, lamlines, mask, tol = tol

if not keyword_set(tol) then tol = 1.
;; tolerance in angstroms. All pixels with values lamlines +/- tol
;; will be masked out

nlines = n_elements(lamlines)

mask = fltarr(n_elements(inlam))
for i = 0, nlines-1 do begin
   
   idx=where(abs(inlam-lamlines[i]) le tol,count)
   if count gt 0 then mask[idx] = 1

endfor
return
end
