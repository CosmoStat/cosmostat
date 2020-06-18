function op_example,vector
	mat=randomn(0,100,200)
	return, mat##vector
end

function op_example_transp,vector
	mat=randomn(0,100,200)
	return, transpose(mat)##vector
end

function check_linear_operator_norm,op_name,transp_op_name,size_input,size_output,nit_max
;Power Method
;see Estimating the matrix p-norm, NJ Higham, Numer. Math. 62:539-555, 1992
;inputs:
;	-op_name: name of the function implementing the operator action on a vector
;	-transp_op_name: name of the function implementing the transpose operator action on a vector
;	-size_input: size of the input vector space 
;	-size_output: size of the output vector space
;output:
;	- lower bound on the matrix 2 norm 
 
init_vec=randomn(seed,1,size_input[0])
vec_2=DBLARR(1,size_input[0])
vec_1=DBLARR(1,size_output[0])

init_vec=init_vec/norm(init_vec,l=2)
for k=0,nit_max-1 do begin
  vec_1 = CALL_FUNCTION(op_name,init_vec)
  vec_1n=norm(reform(vec_1,size_output[0]))
  vec_1p=vec_1/vec_1n
  vec_2=  CALL_FUNCTION(transp_op_name,vec_1p)

  if(norm(reform(vec_2,size_input[0])) lt transpose(vec_2)##init_vec) then break
  init_vec=vec_2/norm(reform(vec_2,size_input[0]))
endfor
if (k eq nit_max) then print,"Max number of iterations reached"
gam=vec_1n
return,gam
end


pro test_linear_operator
	mat=randomn(0,100,200)
	svdc,mat,diag,u_m,v_m
	print,max(diag[*]) ;True max eigenvalue
	print,check_linear_operator_norm("op_example","op_example_transp",100,200,100);Lower bound on the matrix 2 norm

end

