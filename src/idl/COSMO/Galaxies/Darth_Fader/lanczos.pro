function lanczos, x, a

f = a*sin(!pi*x)*sin(!pi*x/a)/(!pi^2*x^2)
idx = where(abs(x) lt 1.e-5,count)

if count gt 0 then f[idx] = 1.

return, f
end
