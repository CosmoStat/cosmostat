
pro mrsp_read_cl,cl4

name = [ 'wmap_comb_tt_powspec_3yr_v2.txt', 'wmap_te_powspec_3yr_v2.txt','wmap_tb_powspec_3yr_v2.txt','wmap_eb_powspec_3yr_v2.txt'] 
cl4 = dblarr(1000,4)

for i =0,3 do begin

tmp = read_ascii(name[i],COMMENT_SYMBOL="#")
hs,tmp
help,tmp.field1
cl4(tmp.field1(0,*),i) = tmp.field1(1,*)
cl4(tmp.field1(0,*),i) = 2*!dpi * cl4(tmp.field1(0,*),i) / (tmp.field1(0,*)*(tmp.field1(0,*)+1))

endfor
end

pro mrsp_read_cl2,cl4


name = [ 'wmap_lcdm_bf_model_yr1_v1.txt']
cl4 = dblarr(1000,4)
tmp = read_ascii(name,COMMENT_SYMBOL="#")


hs,tmp
help,tmp.field1
cl4(tmp.field1(0,*),0) = tmp.field1(1,*)
cl4(tmp.field1(0,*),0) = 2*!dpi * cl4(tmp.field1(0,*),0) / (tmp.field1(0,*)*(tmp.field1(0,*)+1))

cl4(tmp.field1(0,*),1) = tmp.field1(3,*)
cl4(tmp.field1(0,*),1) = 2*!dpi * cl4(tmp.field1(0,*),1) / (tmp.field1(0,*)*(tmp.field1(0,*)+1))


cl4(tmp.field1(0,*),3) = tmp.field1(2,*)
cl4(tmp.field1(0,*),3) = 2*!dpi * cl4(tmp.field1(0,*),3) / (tmp.field1(0,*)*(tmp.field1(0,*)+1))

end


