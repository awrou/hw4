%% Plot Array
File = 'data/temp_array_serial_100.txt';
A = dlmread(File,' ',3,2);
nx=sqrt(length(A));
T=reshape(A,[nx,nx])'

image(T)
image(T,'CDataMapping','scaled')
colorbar
tmp_avg=sum(A)/nx^2