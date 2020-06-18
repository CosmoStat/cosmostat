;+
; NAME:
;        get_stat
;
; PURPOSE:
;		Return statistical information relative to a given data set. The return value is an IDL array of 9 elements.
;			Tab[0] = standard deviation
;			Tab[1] = skewness
;			Tab[2] = Kurtosis
;			Tab[3] = Min
;			Tab[4] = Max 
;			Tab[5] = HC
;			Tab[6] = HC^+
;			Tab[7] = Cumulant of order 5
;			Tab[8] = Cumulant of order 6
;			Tab[9] = Mean
;
;       If the keyword norm is set, the data are first normalized.
;
; CALLING:
;
;		TabStat = get_stat( Data, HCIma=HCIma, TabStatName=TabStatName, norm=norm, qpplot=qpplot, verb=verb, zeromean=zeromean, TabCumulant=TabCumulant ) 
;       
; INPUTS:
;		Data -- IDL array : Input data to analyze
;
; INPUT KEYWORDS:
;		Norm : scalar -- if set, the input data are normalized and centered (i.e. Data = (Data-Mean)/Sigma)
;     qpplot : scalar -- if set, plot the qpplot of the data
;		verb : scalar -- if set, the calculated statists are printed on the screen
;	zeromean : scalar -- if set, the input data are supposed to have a zero mean and are not centered, it is ignored if keyword norm is set
;
; OUTPUT KEYWORDS: 
;     TabStatName -- IDL table of string: TabStatName = [ "Sigma", "Skewness", "Kurtosis", "Min", "Max", "HC1", "HC2", "CUMULANT ORDER 5", "CUMULANT ORDER 6" ]
;     TabCumulant -- IDL double array [0:5]: 6 first cumulants of Data  ( TabCumulant[c] = cumulant of order c+1 )
;
; EXAMPLE:
;       TabStat = get_stat(Data, /verb)
;
; EXTERNAL CALLS:
;
; HISTORY:
;	Written: Jean-Luc Starck, 2005
;	September, 2005 File creation
;-
;=========================================================================

pro HCT, Data, HCIma=HCIma, HCI, HCI2, zeromean=zeromean, verb=verb, nonorm=nonorm

N = N_ELEMENTS(Data)
ninv = double(1.) / double(n)
 
if not keyword_set(zeromean) then m = mean(Data) else m = 0
x = dblarr(N)
x[*] = Data
;s = mad(x)
s = sigma(Data)
if keyword_set(nonorm) then y=x $
else y = (x - m) / s
phi = (double(1.) + double(erf(y/sqrt(2.)))) / 2.
p = 1. - phi
indsort = sort(p)
ps =  p(indsort)
ind = where (ps LT ninv, count)

v = (findgen(N) + 1.) / float(N)
PSN = ps - ps^2.
ind = where ( PSN LT 0, count)
if count GT 0 then PSN[ind] = 0
den = sqrt(PSN)

ind = where (den GT  5.96047e-08, count)
w = ps
w(*) = 0
if count GT 0 then  w(ind) = sqrt( float(N) ) * ( v(ind) - ps(ind)) / den(ind)
HCIma = Data
HCIma[indsort] = w

;HCI = max( w(0:N/2))
HCI = max(abs(w(0:N-1)))
ind = where ( p LT (1./float(N)), count)
a = count
ind = where ( p GT (1- (1./float(N))), count)
b = count
if a GE N/2. then begin
  HCI2 = 0.
;else HCI2 = max(w(a:N/2))
endif else begin
  HCI2 = max(abs(w(a:N-b-1)))
endelse
if keyword_set(verb) then print,  ' HC1 = ', HCI,  ' HC2 = ', HCI2
end

;=========================================================================

function skew, s

n=double(n_elements(s))

ss= skewness(s)*n/(n-1.)*n/(n-2.)

return, ss
end

;=========================================================================

function kurt, s

m=mean(s)
n=double(n_elements(s))

k4=total((s-m)^4)/n
k2=total((s-m)^2)/n

kk=n*n*((n+1)*k4-3*(n-1)*k2*k2)/(n-1)/(n-2)/(n-3)/(k2/n*(n-1))^2
return, kk
end

;=========================================================================


function my_moment, x1, p, opt
x=x1
if  opt EQ  'c' then  x = x - mean(x)  
if  opt EQ  'a' then  x = abs (x) 
m = mean(x^p) 

return, m
end

;=========================================================================

function cumulants, x, p
kp=-1.
;function kp=cumulants(x,p)

; usage:  cumulants (x, p)
; Computes the unbiased estimator of the p-th cumulant up to 6th order of vector
; x using k-statistics, except for orders 5 and 6 where biased estimators are used.
; Fisher R.A., Contributions to mathematical statistics. Wiley and Sons, (1950)

if p GT 6 or  p LT 1 then begin
         print, 'Order of kn must be <=6 and >=1'
	 goto, DONE
end

vs = size(x)
n= double(vs[1])

mc = dblarr(7)
mc(1)=mean(x) 
for k=2,6 do mc(k)=my_moment(x,k,'c')
 
if p EQ 1 then kp=mc(1)
if p EQ 2 then kp=n*mc(2)/(n-1)
if p EQ 3 then kp=n^2*mc(3)/((n-1)*(n-2))
if p EQ 4 then kp=(n^2*((n+1)*mc(4)-3*(n-1)*mc(2)^2))/((n-1)*(n-2)*(n-3))
if p EQ 5 then kp=(n^3*((n+5)*mc(5)-10*(n-1)*mc(2)*mc(3)))/((n-1)*(n-2)*(n-3)*(n-4))
 	; Biased: kp=mc(5,:)-10*mc(2,:)*mc(3,:);
if p EQ 6 then kp=(n^2*((n+1)*(n^2+15*n-4)*mc(6)-15*(n-1)^2*(n+4)*mc(2)*mc(4) $
		-10*(n-1)*(n^2-n+4)*mc(3)^2+30*n*(n-1)*(n-2)*mc(2)^3))/((n-1)*(n-2)*(n-3)*(n-4)*(n-5))
 	; Biased
 	;m(1)=mc(1);
 	;for k=2:6
	;	m(k)=moment(x,k);
	;end
 	;kp=-120*m(1)^6+360*m(1)^4*m(2)-270*m(1)^2*m(2)^2+30*m(2)^3 ...
	;   -120*m(1)^3*m(3)+120*m(1)*m(2)*m(3)-10*m(3)^2+30*m(1)^2*m(4)  ...
	;   -15*m(2)*m(4)-6*m(1)*m(5)+m(6);
 
DONE:
   return, kp
end

;=========================================================================

function get_stat, Data, HCIma=HCIma, TabStatName=TabStatName, norm=norm, qpplot=qpplot, verb=verb, zeromean=zeromean, TabCumulant=TabCumulant
TabStat=0
if N_PARAMS() LT 1 then begin 
        print, 'CALLING SEQUENCE: TabStat = get_stat( Data, Data, HCIma=HCIma, TabStatName=TabStatName, norm=norm, qpplot=qpplot, verb=verb, zeromean=zeromean, TabCumulant=TabCumulant )'
        goto, DONE
        end
 	
N = N_ELEMENTS(Data)
X = dblarr(N)

if keyword_set(norm) then X[*] = (Data - mean(Data)) / sigma(data) $
else X[*] = Data

NStat = 11
TabStatName = ["Sigma", "Skewness", "Kurtosis", "Min", "Max", "HC1", "HC2", "CUMULANT Order 5", "CUMULANT Order 6", "Mean", "MAD"]
TabStat = dblarr(NStat)
TabStat[0] = sigma(X)
TabStat[1] = skew(X)
TabStat[2] = kurt(X)
TabStat[3] = min(X)
TabStat[4] = max(X)
TabStat[9] = mean(X)
TabStat[10] = mad(X)

HCT, x, HCI, HCI2, zeromean=zeromean, HCIma=HCIma
TabStat[5] = HCI
TabStat[6] = HCI2

NbrCumulant=6
TabCumulant = dblarr(NbrCumulant)
if not keyword_set(norm) then x = x / sigma(x)
for c=0,5 do  TabCumulant[c] = cumulants(x, c+1)
TabStat[7] = TabCumulant[4]
TabStat[8] = TabCumulant[5]

if keyword_set(qpplot) then begin
  n = N_ELEMENTS(X)
  s = median(abs(x))/ 0.6745
  y = (x - mean(x)) / s
  noise = randomn(seed, N)
  qpplot, noise, y
end

if keyword_set(verb) then begin
for i=0, NStat-1 do print, '      ', TabStatName[i], ' = ', TabStat[i]
end

DONE:

return, TabStat

end

;=========================================================================


