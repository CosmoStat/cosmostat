function R_to_lstep, R, inverse = inverse

A = 2.*sqrt(2.*alog(2))

if not keyword_set(inverse) then lstep = alog10((A*R + 0.5)/(A*R - 0.5)) else $
   lstep = 0.5*(10^R + 1.)/(A*(10^R - 1.))

return, lstep

end

