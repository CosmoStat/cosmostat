
FUNCTION Spread_Pckt_Tjul, Tjul, Band

; Procedure to spread a 1-d array of science packet Julian times, Tjul, to a
; 2-d array of times evaluated at each observation for a given band.

; Input
;  Tjul		A 1-d array of reduced Julian times corresponding to the time
;		stamps in a given array of science packets.
;  Band		A string containing the band for which to spread the times.
;
; Outputs
;  Tjul_2	A 2-d array of reduced Julian times. Tjul(0,*) is identical to
;		the input array, Tjul(i,*) gives the times at each observation
;		within a packet.

; Set N_obs and do some error checking
Band = STRUPCASE(Band)
CASE Band OF
 'K'  : N_obs = 12L
 'KA' : N_obs = 12L
 'Q'  : N_obs = 15L
 'V'  : N_obs = 20L
 'W'  : N_obs = 30L
 ELSE : MESSAGE,' Invalid band string supplied'
ENDCASE

; Form an array of relative observation times
Tj_obs = (DINDGEN(N_obs)/DOUBLE(N_obs))*1.536d0/86400.d0

; Generate the 2-d array of absolute observation times

N_sci = N_ELEMENTS(Tjul)
Tjul_2 = REPLICATE(1.d0,N_obs)#Tjul + Tj_obs#REPLICATE(1.d0,N_sci)

RETURN, Tjul_2
END
