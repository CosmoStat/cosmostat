function ITrans = DyadDecompSkip(I,J);

[N,T] = size(I);
I00 = I(1:2:N,1:2:T);
Trans00 = cwpt2_interface(I00,'forward','spline','tri',J);
I10 = I(2:2:N,1:2:T);
Trans10 = cwpt2_interface(I10,'forward','spline','tri',J);
I01 = I(1:2:N,2:2:T);
Trans01 = cwpt2_interface(I01,'forward','spline','tri',J);
I11 = I(2:2:N,2:2:T);
Trans11 = cwpt2_interface(I11,'forward','spline','tri',J);
ITrans = zeros(N,T,2*J+1);
ITrans(1:2:N,1:2:T,:)=Trans00(:,:,:);
ITrans(2:2:N,1:2:T,:)=Trans10(:,:,:);
ITrans(1:2:N,2:2:T,:)=Trans01(:,:,:);
ITrans(2:2:N,2:2:T,:)=Trans11(:,:,:);


 
