% Matlab script to compile all functions of the 2D Continuous Wavelet Packet Transform package.
% See 'readme.txt' for more instructions.

if (isempty(who('compile_cwpt2_failed')))      % if not defined...
    %cd cwpt2/
    compile_cwpt2_failed = 1;
end;

% default CFLAGS in Matlab 6: -fPIC -ansi -D_GNU_SOURCE -pthread
% for Matlab 7 one should assure that CFLAGS='-fPIC -ansi -D_GNU_SOURCE -pthread -fexceptions'
% to avoid Matlab closing when getting an exception thrown by mexErrMsgTxt()

mex cwpt2.c dyadic2.c atrous.c btree.c quad.c
mex cwpt2i.c dyadic2.c atrous.c btree.c quad.c
mex c3wpt.c dyadic2.c atrous.c btree.c quad.c
mex c3wpti.c dyadic2.c atrous.c btree.c quad.c
mex cwt.c dyadic2.c atrous.c btree.c quad.c
mex cwti.c dyadic2.c atrous.c btree.c quad.c
mex cwpt2_btree.c btree.c  
mex cwpt2_get_leaf_order.c btree.c

clear compile_cwpt2_failed
cd ..

% 2D Continuous Wavelet Packet Transform package
% (c) 2002-2005 Let It Wave, all rights reserved
