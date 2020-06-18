


;-------------------------------------------------------------------------------
;+
; NAME:
;	jade
;	
; PURPOSE:
;	Apply the Independant Componant Method to a set of vectors using
;       the JADE method.
;
; EXPLANATION:
;	compute ICA
;
; CALLING SEQUENCE:
;	jade, ObservSig, NbSource, DeMixingMat, Process $
;             DispDebug=DispDebug, FileDebug=FileDebug, $
;	      Verbose=Verbose
; INPUTS:
;	ObservSig : mixing of input signal (ObservSig = A # input signals)
;                    ObservSig(i,*) = ith vector
;       NbSource  : number of sources in input signal
;
; OPTIONAL INPUTS:
;	DispDebug : debug trace on screen
;	FileDebug : debug trace on file (Trace.dat in local directory)
;	Verbose	: verbose
;
; KEYWORD PARAMETERS:
;
; OUTPUTS:
;	DeMixingMat : demixing matrix 
;       Process: Reconstructed process = DeMixingMat # ObservSig
;                Process(i, *) = ith Process, with i = 0.. NbSource - 1
;
; DEPENDENCIES:
;
; RESTRICTIONS:
;
; PROCEDURES USED:
;	whitening_signal, estimation_cumulant, init_diagonalisation
;	contrast_optimisation, separating_matrix
;
; EXAMPLE:
;	jade, ObservSig, NbSource, DeMixingMat
;
; MODIFICATION HISTORY:
;
;-------------------------------------------------------------------------------   

;-------------------------------------------------------------------------------
;
; NAME:	
;	whitening_signal
;
; PURPOSE:
;	whiten signal 
;	
; EXPLANATION:
;	remove mean, compute covariance and get the NbSource more significant
;	eigen value to compute the whiten matrix
;
; CALLING SEQUENCE:
;	whitening_signal, NbSensor, NbSource, ObservSig, $
;                         FileDescript, WinNumber, SigLength, $
;                         WhitenObservSig, WhiteningMat, InvWhiteningMat, $
;	                  DispDebug=DispDebug, FileDebug=FileDebug, $
;		          Verbose=Verbose
;
; INPUTS:
;	ObservSig : mixing of input signal (ObservSig = A # input signals)
;       NbSource  : number of source in input signal
;	NbSensor  : number of sensor in observed signal
;	SigLength : length of observed signal
;	FileDescript : file descriptor used if FileDebug is set
;	WinNumber : window number used if DispDebug is set
;
; OPTIONAL INPUTS:
;	DispDebug : debug trace on screen
;	FileDebug : debug trace on file (Trace.dat in local directory)
;	Verbose	  : verbose
;
; KEYWORD PARAMETERS:
;
; OUTPUTS:
;	WhitenObservSig : whiten signal (= whiten matrix # observed signal) 
;	WhiteningMat : whiten matrix
;	InvWhiteningMat : (pseudo)inverse of whiten matrix
;
; DEPENDENCIES:
;
; RESTRICTIONS:
;
; PROCEDURES USED:
;	mean, transpose, eigenql, sqrt, evecs
;
; EXAMPLE:
;
; MODIFICATION HISTORY:
;
;-------------------------------------------------------------------------------
pro whitening_signal, NbSensor, NbSource, ObservSig, $
                      FileDescript, WinNumber, SigLength, $
                      WhitenObservSig, WhiteningMat, InvWhiteningMat, $
	              DispDebug=DispDebug, FileDebug=FileDebug, $
		      Verbose=Verbose


   ;remove mean of observed signal
   ;------------------------------
   if (Verbose) then $
      print, "jade => Removing the mean value"
      
   for curSensor=0,NbSensor-1 do $
      ObservSig[curSensor,*] =   ObservSig[curSensor,*] $
                               - mean(ObservSig[curSensor,*])
			       
   if (FileDebug) then begin
      printf, FileDescript, "Remove mean of observed signal"
      printf, FileDescript, "------------------------------"
      printf, FileDescript, ""
      for curSensor=0,NbSensor-1 do $
         printf, FileDescript, "  Mean value signal[", $
	                       strcompress(curSensor,/remove_all), $
	                       "] : ", mean(ObservSig[curSensor,*])
      printf, FileDescript, ""
   endif
   
   if (DispDebug) then begin
      window, WinNumber, XSize=500, YSize=500, Title='Centered obseved Signal'
      !p.multi=[0,0,NbSensor,0,0]
      for i=0,NbSensor-1 do $
         plot, ObservSig[i,*]
      WinNumber=WinNumber+1
   endif     
   
   ;whitening Observed Signal
   ;-------------------------
   if (Verbose) then $
      print, "jade => Whitening the data"
      
   CovarianceRxx =  (ObservSig # transpose(ObservSig))/SigLength   
                 ;could be compute by correlate(ObservSig,/covariance)
		 ;size = [NbSensor,NbSensor]
   
   ;compute eigen values and eigein vectors
   eigenvalues = eigenql(CovarianceRxx, EIGENVECTORS = evecs, $   
                         RESIDUAL = residual, /ASCENDING)
                 ;we used eigenql function because covariance matrix is
		 ;always symetric.
		 
   ;take the nb source most significant eigen values
   BegMostSignifVal = NbSensor-1 - NbSource+1
   EndMostSignifVal = NbSensor-1
   Scales = sqrt(eigenvalues[BegMostSignifVal:EndMostSignifVal])
                 ;scale vector is sqrt(eigenval) for the nb source
		 ;most significant source. eigneql sort the eigenvalue 
		 ;=> so the m most significant are at the m last one 
		 ;=> we take in [BegMostSignifVal:EndMostSignifVal] in eigenvalues vector
		 
   ;compute whitening matrix 
   WhiteningMat = transpose (evecs[*,BegMostSignifVal:EndMostSignifVal])
   for i=0,NbSource-1 do $
      WhiteningMat[i,*] = WhiteningMat[i,*] * 1./Scales[i]
   ;print, "whiten matrix : ", WhiteningMat
      
   ;compute inverse matrix
   InvWhiteningMat = evecs[*,BegMostSignifVal:EndMostSignifVal]
   for i=0,NbSource-1 do $
      InvWhiteningMat[*,i]=InvWhiteningMat[*,i] * Scales[i]
		 
   ;whitening observed signal
   WhitenObservSig = WhiteningMat # ObservSig
   
   if (DispDebug) then begin
      window, WinNumber, XSize=500, YSize=500, Title='Withen obseved Signal'
      !p.multi=[0,0,NbSource,0,0]
      for i=0,NbSource-1 do $
         plot, WhitenObservSig[i,*]
      WinNumber=WinNumber+1
   endif     
         
   if (FileDebug) then begin
      printf, FileDescript, "whitening Observed Signal"
      printf, FileDescript, "-------------------------"
      printf, FileDescript, ""
      printf, FileDescript, "  Eigen values : size:", $
                            strcompress((size(eigenvalues))(0),/remove_all), " -> [", $
                            strcompress((size(eigenvalues))(1),/remove_all), "]"
      printf, FileDescript, eigenvalues
      printf, FileDescript, "  Eigen vector : size:", $
                            strcompress((size(evecs))(0),/remove_all), " -> [", $
                            strcompress((size(evecs))(1),/remove_all), $
			    ",",  strcompress((size(evecs))(2),/remove_all), "]" 
      printf, FileDescript, evecs
      printf, FileDescript, "  Scales : size:", $
                            strcompress((size(Scales))(0),/remove_all), " -> [", $
			    strcompress((size(Scales))(1),/remove_all), "]"
      printf, FileDescript, Scales		    
      printf, FileDescript, "  Whitening matrix : size:", $
                            strcompress((size(WhiteningMat))(0),/remove_all), " -> [", $
                            strcompress((size(WhiteningMat))(1),/remove_all), $
			    ",",  strcompress((size(WhiteningMat))(2),/remove_all), "]" 
      printf, FileDescript, WhiteningMat
      printf, FileDescript, "  inv Whitening matrix : size:", $
                            strcompress((size(InvWhiteningMat))(0),/remove_all), " -> [", $
                            strcompress((size(InvWhiteningMat))(1),/remove_all), $
			    ",",  strcompress((size(InvWhiteningMat))(2),/remove_all), "]" 
      printf, FileDescript, InvWhiteningMat
      printf, FileDescript, ""
   endif
   
return
end


;-------------------------------------------------------------------------------
;
; NAME:
;	
; PURPOSE:
;
; EXPLANATION:
;
; CALLING SEQUENCE:
;	
; INPUTS:
;
; OPTIONAL INPUTS:
;
; KEYWORD PARAMETERS:
;
; OUTPUTS:
;
; DEPENDENCIES:
;
; RESTRICTIONS:
;
; PROCEDURES USED:
;
; EXAMPLE:
;
; MODIFICATION HISTORY:
;
;-------------------------------------------------------------------------------   
pro estimation_cumulant, NbSensor, NbSource, WhitenObservSig,$
                         FileDescript, WinNumber, SigLength, $
                         NbCumulantMat, CumulantMat, $
	                 DispDebug=DispDebug, FileDebug=FileDebug, $
		         Verbose=Verbose
					
   if (Verbose) then $
      print, "jade => Estimating the cumulant matrices"
   DimSym = (NbSource*(NbSource+1))/2;
   NbCumulantMat = DimSym
   CumulantMat = dblarr(NbSource,NbSource*NbCumulantMat)
   CumulantMat[*,*]=0
   R = dblarr(NbSource,NbSource) 
   R[*,*]=0
   R[indgen(NbSource),indgen(NbSource)]=1.
   Qij = dblarr(NbSource)
   Qij[*]= 0.
   Xim=dblarr(1,NbSource)
   Xim[*,*]=0
   Xjm=dblarr(1,NbSource)
   Xjm[*,*]=0
   Scale = dblarr(NbSource,1) 
   Scale[*,*] = 1./SigLength 
   
   ;a little uncomprehensible trick ......
   X = WhitenObservSig
   IndBeg = 0
   IndEnd = NbSource-1
   for im=0,NbSource-1 do begin
         Xim = X[im,*]
         Qij =   ((Scale # (Xim*Xim)) * X) # transpose(X) $ 
	       - R - 2*(R[*,im]#transpose(R[*,im]))
         CumulantMat[*,IndBeg:IndEnd] = Qij
         IndBeg = IndBeg + NbSource
         IndEnd = IndEnd + NbSource
      for jm=0,im-1 do begin
         Xjm = X[jm,*]
	 Qij =   ((Scale # (Xim*Xjm)) * X) # transpose(X) $
	       - R[*,im]#transpose(R[*,jm]) - R[*,jm]#R[*,im]
	 CumulantMat[*,IndBeg:IndEnd] = sqrt(2.)*Qij
         IndBeg = IndBeg + NbSource
         IndEnd = IndEnd + NbSource  
      endfor
   endfor 
   
   if (FileDebug) then begin
      printf, FileDescript, "Estimation of the cumulant Matrix"
      printf, FileDescript, "---------------------------------"
      printf, FileDescript, ""
      printf, FileDescript, "  DimSym : ", $
                            strcompress(DimSym,/remove_all)
      printf, FileDescript, "  Cumulant matrix number : ", $
                            strcompress(NbCumulantMat,/remove_all)		    
      printf, FileDescript, "  Cumulant matrix : size:", $
                            strcompress((size(CumulantMat))(0),/remove_all), " -> [", $
                            strcompress((size(CumulantMat))(1),/remove_all), $
			    ",",  strcompress((size(CumulantMat))(2),/remove_all), "]" 
      printf, FileDescript, CumulantMat
      printf, FileDescript, ""
   endif   
   			
return
end		
			
			
			
;-------------------------------------------------------------------------------
;
; NAME:
;	
; PURPOSE:
;
; EXPLANATION:
;
; CALLING SEQUENCE:
;	
; INPUTS:
;
; OPTIONAL INPUTS:
;
; KEYWORD PARAMETERS:
;
; OUTPUTS:
;
; DEPENDENCIES:
;
; RESTRICTIONS:
;
; PROCEDURES USED:
;
; EXAMPLE:
;
; MODIFICATION HISTORY:
;
;-------------------------------------------------------------------------------   
pro init_diagonalisation, NbSensor, NbSource, NbCumulantMat, CumulantMat, $
                          FileDescript, WinNumber, SigLength, $
			  eval, evec, $
 			  DispDebug=DispDebug, FileDebug=FileDebug, $
		          Verbose=Verbose 		
			
   if (Verbose) then $
      print, "jade => Initialization of diagonalization"
      			
   ;compute eigen values
   eval = HQR(ELMHES(CumulantMat[*,0:NbSource-1]), /DOUBLE)
   
   ;compute eigin vector
   evec = double(EIGENVEC(CumulantMat[*,0:NbSource-1], eval, RESIDUAL = residual))
   eval = double(eval)
   
   for u=0,NbSource*NbCumulantMat-1,NbSource do $
       CumulantMat[*,u:u+NbSource-1] = CumulantMat[*,u:u+NbSource-1] # evec
   CumulantMat = float (transpose(evec) # CumulantMat)
   
   if (FileDebug) then begin
      printf, FileDescript, "Init joint diagonalization  of the cumulant Matrix"
      printf, FileDescript, "--------------------------------------------------"
      printf, FileDescript, ""
      printf, FileDescript, "  Eigen values : size:", $
                            strcompress((size(eval))(0),/remove_all), " -> [", $
                            strcompress((size(eval))(1),/remove_all), "]"
      printf, FileDescript, eval
      printf, FileDescript, "  Eigen vector : size:", $
                            strcompress((size(evec))(0),/remove_all), " -> [", $
                            strcompress((size(evec))(1),/remove_all), $
			    ",",  strcompress((size(evec))(2),/remove_all), "]" 
      printf, FileDescript, evec
      printf, FileDescript, "  Cumulant matrix : size:", $
                            strcompress((size(CumulantMat))(0),/remove_all), " -> [", $
                            strcompress((size(CumulantMat))(1),/remove_all), $
			    ",",  strcompress((size(CumulantMat))(2),/remove_all), "]" 
      printf, FileDescript, CumulantMat
      printf, FileDescript, ""
   endif	
   
return
end	   



;-------------------------------------------------------------------------------
;
; NAME:
;	
; PURPOSE:
;
; EXPLANATION:
;
; CALLING SEQUENCE:
;	
; INPUTS:
;
; OPTIONAL INPUTS:
;
; KEYWORD PARAMETERS:
;
; OUTPUTS:
;
; DEPENDENCIES:
;
; RESTRICTIONS:
;
; PROCEDURES USED:
;
; EXAMPLE:
;
; MODIFICATION HISTORY:
;
;-------------------------------------------------------------------------------   
pro contrast_optimisation, NbSensor, NbSource, NbCumulantMat, CumulantMat, $
                           evec, FileDescript, WinNumber, SigLength, $
			   DispDebug=DispDebug, FileDebug=FileDebug, $
		           Verbose=Verbose 	
			   
   ;joint diagonalization proper
   if (Verbose) then $
      print, "jade => Contrast optimisation by joint diagonalisation"
   if (FileDebug) then begin
      printf, FileDescript, "Contrast optimisation by joint diagonalisation" 
      printf, FileDescript, "----------------------------------------------"
      printf, FileDescript, ""
   endif    
   
   ;init
   ;----		   
   Seuil = 1./sqrt(SigLength)/100.
   encore = 0
   sweep = 0
   updates = 0
   g = dblarr(2,NbCumulantMat)
   g[*,*] = 0.
   gg = dblarr(2,2)
   gg[*,*] = 0.
   G= dblarr(2,2)
   G[*,*] = 0.
   c = double(0.)
   s = double(0.)
   ton = double(0.)
   toff = double(0.)
   theta = double(0.)
   
   ;optim loop
   ;----------
   while (encore eq 0) do begin
   
      if (Verbose) then $
         print, "jade => Sweep ", strcompress(Sweep,/remove_all)
      if (FileDebug) then $
      	 printf, FileDescript, "  Sweep : ", strcompress(Sweep,/remove_all) 
      Sweep = Sweep + 1
   
      for p=0,NbSource-2 do begin
         for q=p+1,NbSource-1 do begin
	 
	    Ip = indgen((NbSource*NbCumulantMat-1-p)/NbSource+1) * NbSource + p
	    Iq = indgen((NbSource*NbCumulantMat-1-q)/NbSource+1) * NbSource + q
	    
	    ;print, "Ip:", Ip
	    ;print, "Iq:", Iq
	    
	    ;computation of given angle
	    g = [CumulantMat[p,Ip] - CumulantMat[q,Iq], $
	         CumulantMat[p,Iq] + CumulantMat[q,Ip]]
	    gg = g # transpose(g)
	    ton = gg[0,0] - gg[1,1]
	    toff = gg[0,1] + gg[1,0]
	    theta = 0.5 * atan(toff, ton+sqrt(ton*ton+toff*toff))
	    
	    if (FileDebug) then begin
	       printf, FileDescript, "  Ip : ", Ip
	       printf, FileDescript, "  Iq : ", Iq	       
	       printf, FileDescript, "  g : size:", $
                            strcompress((size(g))(0),/remove_all), " -> [", $
                            strcompress((size(g))(1),/remove_all), $
			    ",",  strcompress((size(g))(2),/remove_all), "]" 
               printf, FileDescript, g
	       printf, FileDescript, "  gg : size:", $
                            strcompress((size(gg))(0),/remove_all), " -> [", $
                            strcompress((size(gg))(1),/remove_all), $
			    ",",  strcompress((size(gg))(2),/remove_all), "]" 
               printf, FileDescript, gg
	       printf, FileDescript, "  ton = ", strcompress(ton), $
	                             ", toff = ", strcompress(toff)
	       printf, FileDescript, "  theta = ", strcompress(theta), $
	                             ", Seuil = ", strcompress(Seuil)
	    endif
	    
	    ;givens update
	    if (abs(theta) gt Seuil) then begin
	    
	       Updates = Updates + 1
	       c = cos(theta)
	       s = sin(theta)
	       G = [[c,s],[-s,c]]
	       
	       pair = [[p],[q]]
	       evec [*,pair] = evec[*,pair] # G
	       CumulantMat[pair,*] = transpose(G) # CumulantMat[pair,*]
	       CumulantMat[*,Ip] = c*CumulantMat[*,Ip]+s*CumulantMat[*,Iq]
	       CumulantMat[*,Iq] = -s*CumulantMat[*,Ip]+c*CumulantMat[*,Iq]
		   
	       if (FileDebug) then begin
	          printf, FileDescript, "  Updates : ", $
		                        strcompress(Updates,/remove_all), $
					", c : ", $
					strcompress(c,/remove_all), $
					", s : ", $
					strcompress(s,/remove_all)
		  printf, FileDescript, "  G : size:", $
                            strcompress((size(G))(0),/remove_all), " -> [", $
                            strcompress((size(G))(1),/remove_all), $
			    ",",  strcompress((size(G))(2),/remove_all), "]" 
		  printf, FileDescript, G
		  printf, FileDescript, "  pair : size:", $
                            strcompress((size(pair))(0),/remove_all), " -> [", $
                            strcompress((size(pair))(1),/remove_all), "]" 
		  printf, FileDescript, pair		  	    		
		  printf, FileDescript, "  evec : size:", $
                            strcompress((size(evec))(0),/remove_all), " -> [", $
                            strcompress((size(evec))(1),/remove_all), $
			    ",",  strcompress((size(evec))(2),/remove_all), "]" 
		  printf, FileDescript, evec	
                  printf, FileDescript, "  Cumulant matrix : size:", $
                            strcompress((size(CumulantMat))(0),/remove_all), " -> [", $
                            strcompress((size(CumulantMat))(1),/remove_all), $
			    ",",  strcompress((size(CumulantMat))(2),/remove_all), "]" 
                  printf, FileDescript, CumulantMat		          
	          printf, FileDescript, "  p:", strcompress(p,/remove_all), $
		                        ", q:", strcompress(q,/remove_all), $
			                ", s:", strcompress(s,/remove_all)
					
               endif			       
	 
	    endif else begin 
	       encore = 1
	    endelse
	 endfor
      endfor
   endwhile
   
   if (Verbose) then $
      print, "jade => Total of ", strcompress(Updates,/remove_All), $
             " rotations"
   if (FileDebug) then begin
      printf, FileDescript, "  Total of ", strcompress(Updates,/remove_All), $
                            " rotations"
      printf, FileDescript, "  Eigen vector : size:", $
                            strcompress((size(evec))(0),/remove_all), " -> [", $
                            strcompress((size(evec))(1),/remove_all), $
			    ",",  strcompress((size(evec))(2),/remove_all), "]" 
      printf, FileDescript, evec       
      printf, FileDescript, ""
   endif     		   

return
end			   
			   
			   
			   
;-------------------------------------------------------------------------------
;
; NAME:
;	
; PURPOSE:
;
; EXPLANATION:
;
; CALLING SEQUENCE:
;	
; INPUTS:
;
; OPTIONAL INPUTS:
;
; KEYWORD PARAMETERS:
;
; OUTPUTS:
;
; DEPENDENCIES:
;
; RESTRICTIONS:
;
; PROCEDURES USED:
;
; EXAMPLE:
;
; MODIFICATION HISTORY:
;
;-------------------------------------------------------------------------------   			   
pro separating_matrix, NbSensor, NbSource, evec, WhiteningMat, InvWhiteningMat, $
                       FileDescript, WinNumber, SigLength, $
		       DeMixingMat, $
		       DispDebug=DispDebug, FileDebug=FileDebug, $
		       Verbose=Verbose 		   
			   
   if (Verbose) then $
      print, "jade => Contrast optimisation by joint diagonalisation"		   
			   
   DeMixingMat = transpose(evec) # WhiteningMat
   DeMixingMat = float(DeMixingMat)
   
   if (FileDebug) then begin
      printf, FileDescript, "Separating matrix" 
      printf, FileDescript, "-----------------"
      printf, FileDescript, ""  
      printf, FileDescript, "  Separating matrix : size:", $
                            strcompress((size(DeMixingMat))(0),/remove_all), " -> [", $
                            strcompress((size(DeMixingMat))(1),/remove_all), $
			    ",",  strcompress((size(DeMixingMat))(2),/remove_all), "]"      
      printf, FileDescript, DeMixingMat 
      printf, FileDescript, "" 
   endif
   
   if (Verbose) then $
      print, "jade => Sorting the component"
   A = float(InvWhiteningMat # evec)
   prov = fltarr(NbSource)
   for curSource=0,NbSource-1 do $
      prov[curSource] = total(A[*,curSource] * A[*,curSource])
   SortedInd = sort(prov)
   DeMixingMat = DeMixingMat[SortedInd,*]
   Ind = -1 * indgen(NbSource) + NbSource-1
   DeMixingMat = DeMixingMat[Ind,*]
   
   if (FileDebug) then begin
      printf, FileDescript, "Sorted separating matrix" 
      printf, FileDescript, "------------------------"
      printf, FileDescript, ""  
      printf, FileDescript, "  A matrix : size:", $
                            strcompress((size(A))(0),/remove_all), " -> [", $
                            strcompress((size(A))(1),/remove_all), $
			    ",",  strcompress((size(A))(2),/remove_all), "]"      
      printf, FileDescript, A
      printf, FileDescript, "  vars : size:", $
                            strcompress((size(prov))(0),/remove_all), " -> [", $
                            strcompress((size(prov))(1),/remove_all), $
			    ",",  strcompress((size(prov))(2),/remove_all), "]"      
      printf, FileDescript, prov
      printf, FileDescript, "  key : ", SortedInd
      printf, FileDescript, "  Sorted separating matrix : size:", $
                            strcompress((size(DeMixingMat))(0),/remove_all), " -> [", $
                            strcompress((size(DeMixingMat))(1),/remove_all), $
			    ",",  strcompress((size(DeMixingMat))(2),/remove_all), "]"      
      printf, FileDescript, DeMixingMat
      printf, FileDescript, "" 
   endif  
   
   
   ;Sign are fixed by forcing the firts column of B to have 
   ;non negative entries
   if (Verbose) then $
      print, "jade => Fixing the signs"
   inter = DeMixingMat[*,1]
   IndPos = where (inter ge 0, CountPos)
   IndNeg = where (inter lt 0, CountNeg)
   prov = fltarr(NbSource,NbSource)
   prov[*,*] = 0.
   if (CountPos gt 0) then prov[IndPos,IndPos] =  1
   if (CountNeg gt 0) then prov[IndNeg,IndNeg] = -1
   DeMixingMat = prov # DeMixingMat
   
   if (FileDebug) then begin
      printf, FileDescript, "Fixed sign sorted separating matrix" 
      printf, FileDescript, "-----------------------------------"
      printf, FileDescript, ""  
      printf, FileDescript, "  Fixed sign sorted separating matrix : size:", $
                            strcompress((size(DeMixingMat))(0),/remove_all), " -> [", $
                            strcompress((size(DeMixingMat))(1),/remove_all), $
			    ",",  strcompress((size(DeMixingMat))(2),/remove_all), "]"      
      printf, FileDescript, DeMixingMat
      printf, FileDescript, "" 
   endif    
   
return
end
   			   
;============================================
;============================================
;============================================
		   
pro jtest1, ps=ps

   ;init
   ;----
   NbSource = 3         ; number of sources
   NbSensor = 6        ; number of sensors
   SigLength = 200      ; sample size
   NiseLeveldB = -20    ; noise level in dB
   seed = 37
   
   FirstTime = 1
   Seed = 3764524
   
   begin
   
         freq1 = 0.013   ; freq source 1
         freq2 = 0.02    ; freq source 2
	 FirstTime = 0
  
      print, "freq1 = " , freq1, "  freq2 = " , freq2
      t = indgen(SigLength)+1
   
      ;source signal
      ;-------------
      s1 =       double(cos (2*!pi*freq1*t))
      s2 =       double(cos (2*!pi*freq2*t))
      IndPos  = where(s2 gt 0, countPos)
      IndNull = where(s2 eq 0, countNull)
      IndNeg  = where(s2 lt 0, countNeg)
      if (countPos gt 0) then s2[IndPos]=1.01
      if (countNull gt 0) then s2[IndNull]=0
      if (countNeg gt 0) then s2[IndNeg]=-1.01
      s3 =       randomn(Seed,SigLength)
      SigVector = [ transpose([s1]),transpose([s2]),transpose([s3])]
   
      ;plot source
      ;-----------
      if not keyword_set(ps) then $ 
                    window, 1, XSize=500, YSize=500, Title='Source signal' $
      else setps, /portrait, filename='fig_ica_source.ps'
      !p.multi=[0,0,NbSource,0,0]
      for i=0,NbSource-1 do $
 	 plot, SigVector[i,*] , yticks=3,   $
         ticklen=0.05,yminor=1,xminor=2, xcharsize=1.5, ycharsize = 1.5, $
	 thick=1, xticklen=0.1, xthick=2, ythick=2
	 
       if keyword_set(ps) then endps
     
      ;mixing and noising source
      ;-------------------------
      MixingMat = randomn(Seed, NbSensor, NbSource)  ;rand mixing matrice
      print, "Mixing matrice :"
      print, MixingMat  
   
      NoiseAmp = 10^(NiseLeveldB/20)
      ObservSig = MixingMat # SigVector ; + NoiseAmp*randomn(Seed, NbSensor, SigLength)
   
      ;plot obsevevd signal
      ;--------------------
      if not keyword_set(ps) then $ 
          window, 2, XSize=500, YSize=500, Title='Observed signal' $
      else setps, /portrait, filename='fig_ica_obs_signal.ps'
      !p.multi=[0,0,NbSensor,0,0]
      for i=0,NbSensor-1 do $
 	 plot, ObservSig[i,*] , yticks=3,   $
         ticklen=0.05,yminor=1,xminor=2, xcharsize=1.5, ycharsize = 1.5, $
	 thick=1, xticklen=0.1, xthick=2, ythick=2
	 
      if keyword_set(ps) then endps
      
      ;separation
      ;----------
      jade, ObservSig, NbSource, DeMixingMat
   
      ;identified signal
      ;-----------------
      IdentSig = DeMixingMat # ObservSig
   
      ;plot demixing signal
      ;--------------------
      if not keyword_set(ps) then $ 
          window, 3, XSize=500, YSize=500, Title='Demixing signal' $
      else setps, /portrait, filename='fig_ica_demixing.ps'
      !p.multi=[0,0,NbSource,0,0]
      for i=0,NbSource-1 do $
 	 plot, IdentSig[i,*] , yticks=4,   $
         ticklen=0.05,yminor=2,xminor=2, xcharsize=1.5, ycharsize = 1.5, $
	 thick=1, xticklen=0.1, xthick=2, ythick=2
	 
      if keyword_set(ps) then endps
     
      ;plot demixing signal
      ;--------------------
      ;window, 4, XSize=500, YSize=500, Title='Diff Demixing signal'
      ;!p.multi=[0,0,NbSource,0,0]
      ;for i=0,NbSource-1 do $
      ;   plot, IdentSig[i,*]- SigVector[i,*]  
   
      ;performance
      ;-----------
      print, "Globa system:"
      print, "If this matrix is close to a product Permutation*Diagonal,"
      print, "then separation was successful." 
      Perf = (DeMixingMat#MixingMat)^2.
      print, float(perf)
     
      
   end 
   
return
end

			   
;============================================			   
pro jtest, ps=ps

   ;init
   ;----
   NbSource = 3         ; number of sources
   NbSensor = 6        ; number of sensors
   SigLength = 200      ; sample size
   NiseLeveldB = -20    ; noise level in dB
   seed = 37
   
   FirstTime = 1
   Seed = 3764524
   
   repeat begin
   
      if (FirstTime eq 1) then begin
         freq1 = 0.013   ; freq source 1
         freq2 = 0.02    ; freq source 2
	 FirstTime = 0
      endif else begin
         freq1 = 0.1 * randomu(Seed)
         freq2 = 0.1 * randomu(Seed)
      endelse
      print, "freq1 = " , freq1, "  freq2 = " , freq2
      t = indgen(SigLength)+1
   
      ;source signal
      ;-------------
      s1 =       double(cos (2*!pi*freq1*t))
      s2 =       double(cos (2*!pi*freq2*t))
      IndPos  = where(s2 gt 0, countPos)
      IndNull = where(s2 eq 0, countNull)
      IndNeg  = where(s2 lt 0, countNeg)
      if (countPos gt 0) then s2[IndPos]=1.01
      if (countNull gt 0) then s2[IndNull]=0
      if (countNeg gt 0) then s2[IndNeg]=-1.01
      s3 =       randomn(Seed,SigLength)
      SigVector = [ transpose([s1]),transpose([s2]),transpose([s3])]
   
      ;plot source
      ;-----------
      if not keyword_set(ps) then $ 
                    window, 1, XSize=500, YSize=500, Title='Source signal' $
      else setps, /portrait, filename='fig_ica_source.ps'
      !p.multi=[0,0,NbSource,0,0]
      for i=0,NbSource-1 do $
 	 plot, SigVector[i,*] , yticks=3,   $
         ticklen=0.05,yminor=1,xminor=2, xcharsize=1.5, ycharsize = 1.5, $
	 thick=1, xticklen=0.1, xthick=2, ythick=2
	 
       if keyword_set(ps) then endps
     
      ;mixing and noising source
      ;-------------------------
      MixingMat = randomn(Seed, NbSensor, NbSource)  ;rand mixing matrice
      print, "Mixing matrice :"
      print, MixingMat  
   
      NoiseAmp = 10^(NiseLeveldB/20)
      ObservSig = MixingMat # SigVector ; + NoiseAmp*randomn(Seed, NbSensor, SigLength)
   
      ;plot obsevevd signal
      ;--------------------
      if not keyword_set(ps) then $ 
          window, 2, XSize=500, YSize=500, Title='Observed signal' $
      else setps, /portrait, filename='fig_ica_obs_signal.ps'
      !p.multi=[0,0,NbSensor,0,0]
      for i=0,NbSensor-1 do $
 	 plot, ObservSig[i,*] , yticks=3,   $
         ticklen=0.05,yminor=1,xminor=2, xcharsize=1.5, ycharsize = 1.5, $
	 thick=1, xticklen=0.1, xthick=2, ythick=2
	 
      if keyword_set(ps) then endps
      
      ;separation
      ;----------
      jade, ObservSig, NbSource, DeMixingMat
   
      ;identified signal
      ;-----------------
      IdentSig = DeMixingMat # ObservSig
   
      ;plot demixing signal
      ;--------------------
      if not keyword_set(ps) then $ 
          window, 3, XSize=500, YSize=500, Title='Demixing signal' $
      else setps, /portrait, filename='fig_ica_source.ps'
      !p.multi=[0,0,NbSource,0,0]
      for i=0,NbSource-1 do $
 	 plot, IdentSig[i,*] , yticks=3,   $
         ticklen=0.05,yminor=2,xminor=2, xcharsize=1.5, ycharsize = 1.5, $
	 thick=1, xticklen=0.1, xthick=2, ythick=2
	 
      if keyword_set(ps) then endps
     
      ;plot demixing signal
      ;--------------------
      ;window, 4, XSize=500, YSize=500, Title='Diff Demixing signal'
      ;!p.multi=[0,0,NbSource,0,0]
      ;for i=0,NbSource-1 do $
      ;   plot, IdentSig[i,*]- SigVector[i,*]  
   
      ;performance
      ;-----------
      print, "Globa system:"
      print, "If this matrix is close to a product Permutation*Diagonal,"
      print, "then separation was successful." 
      Perf = (DeMixingMat#MixingMat)^2.
      print, float(perf)
     
      print, "Hit q for exit"
      print, "Hit any other key for another experiment"
      StopLoop = get_kbrd(1)
      
   endrep until StopLoop eq 'q'
   
return
end
	   
;===================================================================

pro jade, ObservSig, NbSource, DeMixingMat, Process, $
          DispDebug=DispDebug, FileDebug=FileDebug, $
	  Verbose=Verbose

   ;init param in
   ;-------------
   if not keyword_set(DispDebug) then DispDebug=0
   if not keyword_set(FileDebug) then FileDebug=0
   if not keyword_set(Verbose) then Verbose=0    ;always ok...
   WinNumber = 10
   
    ;init
   ;----
   NbSensor = (size(ObservSig))(1)
   SigLength = (size(ObservSig))(2)
   if (NbSource gt NbSensor) then NbSource=NbSensor  
   
   ;write if debug
   ;--------------
   if (FileDebug) then begin
      FileDescript = 1
      Openw, FileDescript, 'Trace.dat'
      printf, FileDescript, "INPUT DATA"
      printf, FileDescript, "----------"
      printf, FileDescript, ""
      printf, FileDescript, "  Nb Source : ", $
                            strcompress(NbSource,/remove_all)
      printf, FileDescript, "  Nb Sensor : ", $
                            strcompress(NbSensor,/remove_all)
      printf, FileDescript, "  Length signal : ", $
                            strcompress(NbSource,/remove_all)
      for i=0,(size(ObservSig))(1)-1 do $
         printf, FileDescript, "  Signal[", strcompress(i,/remove_all), $
	                       "] : ", transpose(ObservSig[i,*])
      printf, FileDescript, ""
   endif
   

   ;whitening signal
   ;----------------
   whitening_signal, NbSensor, NbSource, ObservSig,$
                     FileDescript, WinNumber, SigLength, $
                     WhitenObservSig, WhiteningMat, InvWhiteningMat, $
		     DispDebug=DispDebug, FileDebug=FileDebug, $
		     Verbose=Verbose

   
   ;Estimation of the cumulant Matrix
   ;---------------------------------
   estimation_cumulant, NbSensor, NbSource, WhitenObservSig,$
                        FileDescript, WinNumber, SigLength, $
                        NbCumulantMat, CumulantMat, $
			DispDebug=DispDebug, FileDebug=FileDebug, $
		        Verbose=Verbose
			

   ;joint diagonalisation of the cumulant matrix
   ;--------------------------------------------
   
   ;init by diagonalizing a single cumulant 
   ;---------------------------------------
   init_diagonalisation, NbSensor, NbSource, NbCumulantMat, CumulantMat, $
                         FileDescript, WinNumber, SigLength, $
			 eval, evec, $
 			 DispDebug=DispDebug, FileDebug=FileDebug, $
		         Verbose=Verbose 
			  
			  
   ;Contrast optimisation by joint diagonalisation
   ;----------------------------------------------
   contrast_optimisation, NbSensor, NbSource, NbCumulantMat, CumulantMat, $
                          evec, FileDescript, WinNumber, SigLength, $
			  DispDebug=DispDebug, FileDebug=FileDebug, $
		          Verbose=Verbose 
			  
         

   
   ;Separating matrix
   ;-----------------
   separating_matrix, NbSensor, NbSource, evec, WhiteningMat, InvWhiteningMat, $
                      FileDescript, WinNumber, SigLength, $
		      DeMixingMat, $
		      DispDebug=DispDebug, FileDebug=FileDebug, $
		      Verbose=Verbose 
   
   
   ;close file if debug
   ;-------------------
   if (FileDebug) then close, FileDescript

   process = DeMixingMat # ObservSig
   
return
end


