; PRO odeint,ystart,x1,x2,eps,h1,hmin,nok,nbad,derivs,xp,yp,$ 
;    dxsav=dxsav,hwant=hwant,verbose=verbose 
; 
; Jan 05 - downloaded by A. Refregier from  
;          http://uw.physics.wisc.edu/~forest/Physics311/homework.htm 
;          and comments added 
; Jan 05 - Modified by AR to simplify output of intermediate steps
; PURPOSE: Integration of a system of first order ordinary 
;          differential equation dy_i/dx(x)=f_i(x), i=0,nvar using the 
;          Runge-Kunta method with adaptative step size as described 
;          in Numerical Recipies, section 15.2 
;          This IDL routine was downloaded from the above cite and 
;          was presumably translated from the Numerical Recipies 
;          Fortran or C code 
; INPUT: ystart: vector of starting values 
;        x1,x2: integration range (from x1 to x2) 
;        eps: integration accuracy 
;        h1: guess for first step size 
;        hmin: minimum allowed step size (can be zero) 
;        nok, nbad: number of good and bad (but retried and fixed) steps taken 
;        derivs: user supplied subroutine to compute the derivatives 
;                f_i(x,y_j) with call given by derivs(x,y) 
; KEYWORDS: dxsav: minimum step size to store intermediate value 
;           hwant: desired (i.e. maximum) step size. This is useful
;                  to produce a tabulated solution of intermediate values 
;           verbose: print comments 
; OUTPUT: ytstart contains the value of yi at the end of the integration
;         xp, yp: contains intermediate values if the integration (useful
;                 in conjunction to hwant to produce a tabulated solution)
 
function sgn2a,x,y 
 
if y ne 0.0 then s = y/abs(y) else s = 1. 
return,s*x 
 
end 
 
; $Id: odeint.pro,v 1.2 2000/08/07 20:41:23 aake Exp $ 
 
pro rkck,y,dydx,x,h,yout,yerr,derivs 
 
;      INTEGER n,NMAX 
;      REAL h,x,dydx(n),y(n),yerr(n),yout(n) 
;      EXTERNAL derivs 
;      PARAMETER (NMAX=50) 
;CU    USES derivs 
;      INTEGER i 
;      REAL ak2(NMAX),ak3(NMAX),ak4(NMAX),ak5(NMAX),ak6(NMAX), 
;     *ytemp(NMAX) 
 
A2 =.2  &  A3 =.3      &  A4 =.6      &  A5 =1.  &  A6 =.875 
B21=.2  &  B31=3./40.  &  B32=9./40.  &  B41=.3  &  B42=-.9 
B43=1.2 &  B51=-11./54.  &  B52=2.5  &  B53=-70./27.  &  B54=35./27. 
B61=1631./55296.  &  B62=175./512.  &  B63=575./13824. 
B64=44275./110592.  &  B65=253./4096. 
C1=37./378.  &  C3=250./621.  &  C4=125./594.  &  C6=512./1771. 
DC1=C1-2825./27648.  &  DC3=C3-18575./48384. 
DC4=C4-13525./55296.  &  DC5=-277./14336.  &  DC6=C6-.25 
 
ytemp = y+B21*h*dydx 
ak2   = call_function(derivs,x+A2*h,ytemp) 
ytemp = y+h*(B31*dydx+B32*ak2) 
ak3   = call_function(derivs,x+A3*h,ytemp) 
ytemp = y+h*(B41*dydx+B42*ak2+B43*ak3) 
ak4   = call_function(derivs,x+A4*h,ytemp) 
ytemp = y+h*(B51*dydx+B52*ak2+B53*ak3+B54*ak4) 
ak5   = call_function(derivs,x+A5*h,ytemp) 
ytemp = y+h*(B61*dydx+B62*ak2+B63*ak3+B64*ak4+B65*ak5) 
ak6   = call_function(derivs,x+A6*h,ytemp) 
yout  = y+h*(C1*dydx+C3*ak3+C4*ak4+C6*ak6) 
yerr  = h*(DC1*dydx+DC3*ak3+DC4*ak4+DC5*ak5+DC6*ak6) 
 
end 
;------------------------------------------------------------- 
 
 
pro rkqs,y,dydx,x,htry,eps,yscal,hdid,hnext,derivs 
;      INTEGER n,NMAX 
;      REAL eps,hdid,hnext,htry,x,dydx(n),y(n),yscal(n) 
;      EXTERNAL derivs 
;      PARAMETER (NMAX=50) 
;CU    USES derivs,rkck 
;      INTEGER i 
;      REAL errmax,h,htemp,xnew,yerr(NMAX),ytemp(NMAX),SAFETY,PGROW, 
;     *PSHRNK,ERRCON 
;      PARAMETER (SAFETY=0.9,PGROW=-.2,PSHRNK=-.25,ERRCON=1.89e-4) 
SAFETY=0.9 
PGROW =-.2 
PSHRNK=-.25 
ERRCON=1.89e-4 
 
h=htry 
errflag = 1 
while (errflag) do begin 
  rkck,y,dydx,x,h,ytemp,yerr,derivs   
  errmax=max(abs(yerr/yscal)) 
  errmax=errmax/eps 
;  print,'rkqs: x,h,yerr,yscale,errmax:',x,h,yerr(0),yerr(1),$
;    yscal(0),yscal(1),errmax
  if (errmax gt 1.) then begin  ; error too large; reduce step size
    htemp=SAFETY*h*(errmax^PSHRNK) 
    h=sgn2a(max([abs(htemp),0.1*abs(h)]),h) ; no more than a factor of 10
;    print,'rkqs: reducing step size to h:',h
    xnew=x+h 
    if (xnew eq x) then message,'stepsize underflow in rkqs'
  endif else $                  ; step succeded
    errflag = 0 
end 
if (errmax gt ERRCON) then  $   ; step succeded; compute size of next step
  hnext=SAFETY*h*(errmax^PGROW)  $ 
else  $ 
  hnext=5.*h                    ; not more than a factor of 5 increase
hdid=h 
x=x+h 
y=ytemp 
;print,'rkqs: hdid, hnext:',hdid,hnext
 
end 
 
;--------------------------------------------------------------- 
pro odeint,ystart,x1,x2,eps,h1,hmin,nok,nbad,derivs,xp,yp,$ 
    dxsav=dxsav,hwant=hwant,verbose=verbose 
;      INTEGER nbad,nok,nvar,KMAXX,MAXSTP,NMAX 
;      REAL eps,h1,hmin,x1,x2,ystart(nvar),TINY 
;      EXTERNAL derivs,rkqs 
;      PARAMETER (MAXSTP=10000,NMAX=50,KMAXX=200,TINY=1.e-30) 
;      INTEGER i,kmax,count,nstp 
;      REAL dxsav,h,hdid,hnext,x,xsav,dydx(NMAX),xp(KMAXX),y(NMAX), 
;     *yp(NMAX,KMAXX),yscal(NMAX) 
;      COMMON /path/ kmax,count,dxsav,xp,yp 
if n_elements(ystart) eq 0 then begin 
  print,'odeint,ystart,x1,x2,eps,h1,hmin,nok,nbad,derivs,' 
  print,'       dxsav=dxsav,xp=xp,yp=yp,count=count,hwant=hwant,verbose=verbose' 
  return 
endif 

; declarations 
MAXSTP=10000 
TINY=1.e-30 
n_y=n_elements(ystart)
if n_elements(verbose) eq 0 then $ 
  verbose=MAXSTP+1 $ 
else verbose=max([10,verbose]) 
if (n_elements(dxsav) eq 0) then dxsav=hmin 
if (n_elements(hwant) eq 0) then hwant=abs(2.*(x2-x1)) 
hwant = max([hmin,hwant]) 

; initialise integration variables 
x=x1 
h = 1.0 
h=sgn2a(h1,x2-x1) 
signh=sgn2a(1.,x2-x1) 
nok=0 
nbad=0 
count=0
y=ystart 
xsav=x-2.*dxsav

;; initialise intermediate value arrays
;xp=0.            
;yp=fltarr(n_y)
;xsav=x-2.*dxsav

;kmax=0 
;if (n_elements(xp) ne 0) then begin 
;    if ((size(xp))(0) eq 1) then begin 
;        kmax=(size(xp))(1) 
;        yp=fltarr((size(ystart))(1),kmax) 
;    endif else begin 
;        print,'odeint: xp should be 1-d array' 
;        return 
;    endelse 
;endif 
;if (kmax gt 0) then xsav=x-2.*dxsav 

; integrate 
for nstp=1,MAXSTP do begin 
    ; compute derivatives and compute scaling to monitor accuracy
    dydx = call_function(derivs,x,y) 
    yscal=abs(y)+abs(h*dydx)+TINY 
;    print,'odeint:y,dydx,h,yscal:',y(0),y(1),dydx(0),dydx(1),h,yscal(0),yscal(1)    

    ; store intermediate values
    if (abs(x-xsav) gt abs(dxsav)) then begin 
      count=count+1 
;      xp(count)=x 
;      yp(*,count)=y 
;      xsav=x 
      if count eq 1 then begin
        xp=x
        yp=y
      endif else begin
        xp=[xp,x]
        yp=reform([reform(yp,2*(count-1)),y],2,count)
      endelse
    endif 
 
    ; avoid overshooting and compute next step
    if ((x+h-x2)*(x+h-x1) gt 0.) then h=x2-x 
    rkqs,y,dydx,x,h,eps,yscal,hdid,hnext,derivs 
    if (hdid eq h) then $ 
      nok=nok+1 $ 
    else $ 
      nbad=nbad+1 
    
    ; are we done?
    if ((x-x2)*(x2-x1) ge 0.) then begin 
        ystart=y 
        count=count+1 
;         xp(count)=x 
;         yp(*,count)=y 
        if count eq 1 then begin
          xp=x
          yp=y
        endif else begin
          xp=[xp,x]
          yp=reform([reform(yp,2*(count-1)),y],2,count)
        endelse
        return     ; normal exit
    endif 

    ; check if step size smaller than minimum
    if (abs(hnext) lt hmin) then $ 
      stop,'stepsize smaller than minimum in odeint' 
;    h=hnext 
    h=signh*min(abs([hnext,hwant])) 
    if (nstp mod verbose eq 0) then $ 
      print,'odeint: step ',nstp,' stepsize ',h 
end 

; check if max number of steps was reached
stop,'too many steps in odeint' 

return 
end



