map = randomn(seed, 256, 512)
W = BWT01_DIRECT( map, 3, M)
WI =  BWT01_INVERSE( W, 3)
IF mean( map- WI) le 1.e-9 THEN print, 'OK   ' 
