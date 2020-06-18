;
; Compute the inversed Chi-Square CDF using the GSL function
;
; x: the probability
; n: the number degree of freedom
;
function chi2inv,x,n

 command = './chi2inv '+string(double(x))+' '+string(double(n))+'>'+'toto.dat'
 spawn, command
 
; Open file 
 openr,file,'toto.dat',/GET_LUN

; Read Value
 res=double(0.0)
 readf,file,res

; Close file
 free_lun,file
 spawn,'rm -f toto.dat'

 return,res[0]

end
