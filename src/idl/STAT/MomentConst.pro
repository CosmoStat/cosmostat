function MomentConst,Xin,K1,K2,K3,K4,Nmax

Xout = Xin;
p3 = 0.001;
p4 = 0.001;


for ll=0,Nmax-1 do begin

	;---- Kurtosis Const
	
	eK2 = sigma(Xout)
	eK1 = mean(Xout)
	Xout = K2*(Xout - eK1)/eK2+ K1
	
	eK4 = kurtosis(Xout)
	dJ4 = 8./K2^4.*(Xout - K1)^3.*(eK4 - K4)


	Xout = Xout - p4*dJ4

	;---- Skewness Const
	
	eK2 = sigma(Xout)
	eK1 = mean(Xout)
	Xout = K2*(Xout - eK1)/eK2 + K1
	
	eK3 = skewness(Xout)
	dJ3 = 6./K2^3.*(Xout - K1)^2.*(eK3 - K3)

	
	Xout = Xout - p3*dJ3

	

endfor

return,Xout

end