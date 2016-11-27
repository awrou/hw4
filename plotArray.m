%% Plot Array
File = 'data/temparray_serial_512.txt';
A = dlmread(File,' ',3,2);
nx=ceil(sqrt(length(A)));
T=reshape(A,[nx,nx])'

image(T)
image(T,'CDataMapping','scaled')
colorbar
tmp_avg=sum(A)/nx^2;
title(['Serial nx=' num2str(nx) ' t\_avg=' num2str(tmp_avg)]);
