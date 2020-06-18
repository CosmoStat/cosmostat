;
; Compute the Psi polygamma function using the GSL function
;
; x: the degree of the derivative
; n: the value
;
function psi,x,n

 command = './psi '+string(double(x))+' '+string(double(n))+'>'+'toto.dat'
 spawn, command
 
; Open file 
 openr,file,'toto.dat',/GET_LUN

; Read Value
 res=double(0.0)
 readf,file,res

; Close file and remove
 free_lun,file
 spawn,'rm -f toto.dat'

 return,res[0]

end
