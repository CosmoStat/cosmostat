function repmat3D,M,x,y,z

	A=M
	MM=M

	for i=1,x-1 do begin
		A=[A,M]
		MM=[MM,M]
	endfor

	for i=1,y-1 do begin
		A=[[A],[MM]]
	endfor
	MM=A
	for i=1,z-1 do begin
		A=[[[A]],[[MM]]]
	endfor

	return,A
	
end
