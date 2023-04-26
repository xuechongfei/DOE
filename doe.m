%% 理论输出图
image = im2double(rgb2gray(imread('matlab.jpg')));
image0 = imresize(image,[1024,1024],'bilinear');
idea_image0 = 1-image0;
figure(1)
subplot(1,3,1)
imshow(idea_image0);%显示理论输出光强分布
set(gca,'fontname','times new roman','fontsize',16);    %设置图形对象属性
title('输出面理论光强分布','fontname','华文中宋','fontsize',16);


%% 入射高斯光束参数设置
Nxy = 1024;          %x,y坐标采样点
lambda = 632.8e-9;    %波长为 632nm
k = 2*pi/lambda;    %波数
w0 = 6e-3;        %光斑尺寸 6mm
z = 0;
Z_R = pi*w0^2/lambda;           %瑞利长度
w_z = w0*sqrt(1+(z/Z_R)^2);     %光束在z位置的半径
%坐标设置
[x,y] = meshgrid(linspace(-10e-3,10e-3,Nxy));     %近场直角坐标
E = w0/w_z*exp(-(x.^2+y.^2)/w0^2)*exp(-1i*(k*z));%z=0束腰处入射光场分布
I = (abs(E)).^2;%入射光束腰面光强
phase = angle(E);%入射光束腰面相位
%% 经过DOE元件附加随机初相位
add_phase = (2*rand(Nxy,Nxy) - 1)*pi;
estimated_image = E.*exp(1i*add_phase); %模拟初次经过DOE后高斯光束

%% 迭代计算
tic
idea_image = sqrt(idea_image0);%理论振幅图
for m = 1:100
    F_image = fft2(estimated_image);%第k次迭代时，输出光场
    absF_image = abs(F_image)/max(max(abs(F_image)));
    sim = corrcoef(idea_image,absF_image);%振幅的相关系数
    if sim>=0.98
        break
    end  
    F_image_phase = angle(F_image);%保留相位
    F_estimated_image = idea_image.*exp(1i*F_image_phase);%对输出光场振幅和理论输出光场振幅进行交换
    IF_image = ifft2(F_estimated_image);%反傅里叶变换后，计算的k+1次输入光场
    IF_image_phase = angle(IF_image);%保留相位
    estimated_image = E.*exp(1i*IF_image_phase);%实际第k+1次输入光场
end
toc

subplot(1,3,2)
outimage = abs(F_image);
outimage = outimage/max(max(outimage));
imshow(outimage);
set(gca,'fontname','times new roman','fontsize',16);    %设置图形对象属性
title('输出面振幅分布','fontname','华文中宋','fontsize',16);

subplot(1,3,3)
imshow(IF_image_phase)
set(gca,'fontname','times new roman','fontsize',16);    %设置图形对象属性
title('DOE相位分布','fontname','华文中宋','fontsize',16);


idea_image6 = imresize(image,[128,128]);
idea_image6 = sqrt(idea_image6);
F_image6 = imresize(F_image,[128,128]);%kron(IF_image_phase - rand(64,64) * 0.01, ones(2));
F_image_phase6 = angle(F_image6);
F_estimated_image6 = idea_image6.*exp(1i*F_image_phase6);%对输出光场振幅和理论输出光场振幅进行交换
IF_image6 = ifft2(F_estimated_image6);%反傅里叶变换后，计算的k+1次输入光场

IF_image_phase6 = angle(IF_image6);
[x,y] = meshgrid(linspace(-10e-3,10e-3,128));     %近场直角坐标
E6 = w0/w_z*exp(-(x.^2+y.^2)/w0^2)*exp(-1i*(k*z));%z=0束腰处入射光场分布
estimated_image6 = E6.*exp(1i*IF_image_phase6);%实际第k+1次输入光场
F_image6 = fft2(estimated_image6);
figure(2)
subplot(1,3,1)
imshow(idea_image6)

subplot(1,3,2)
outimage6 = abs(F_image6);
outimag6e = outimage6/max(max(outimage6));
imshow(outimag6e);%显示理论输出光强分布

subplot(1,3,3)
imshow(IF_image_phase6)


