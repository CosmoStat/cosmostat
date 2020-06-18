function mad, data
  n = n_elements(data)
  d = fltarr(n)
  d[*] = abs(data-median(data))
  ind = sort(d)
  d = d[ind]
  return, d[n/2] / 0.6747
end
