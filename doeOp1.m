%% 理论输出图
image = im2double(255 - rgb2gray(imread('1.jpg')));
image0 = imresize(image,[512,512],'bilinear');
idea_image0 = image0;
figure(2)
subplot(1,3,1)
imshow(idea_image0);%显示理论输出光强分布
set(gca,'fontname','times new roman','fontsize',16);    %设置图形对象属性
title('输出面理论光强分布','fontname','华文中宋','fontsize',16);


%% 入射高斯光束参数设置
Nxy = 512;          %x,y坐标采样点
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
simlist = zeros(1,200);
idea_image_bw=1 - im2bw(idea_image0, 0);

P1 = idea_image/sum(sum(idea_image));
seta = 2;
for m = 1:200
    F_image = fft2(estimated_image);%第k次迭代时，输出光场
    absF_image = abs(F_image)*2/Nxy;
    sim = corrcoef(idea_image,absF_image);%振幅的相关系数
    simlist(m) = sim(1,2);
    if sim(1,2)>=0.99
        break
    end  
    


    P2 = absF_image/sum(sum(absF_image));
    
    F_image_phase = angle(F_image);%保留相位
    if m == 1
        F_image_phase1 = F_image_phase;
    end 
    
    leta = ((P1./P2).^seta);
    F_estimated_image = (idea_image + absF_image.*leta).*exp(1i*F_image_phase);%对输出光场振幅和理论输出光场振幅进行交换
    IF_image = ifft2(F_estimated_image);%反傅里叶变换后，计算的k+1次输入光场
    IF_image_phase = angle(IF_image);%保留相位
    estimated_image = E.*exp(1i*IF_image_phase);%实际第k+1次输入光场
end
toc

subplot(1,3,2)
outimage = abs(F_image);
outimage = outimage*2/Nxy;
% outimage = (outimage.^2);
imshow(outimage/max(max(outimage)));
set(gca,'fontname','times new roman','fontsize',16);    %设置图形对象属性
title('输出面振幅分布','fontname','华文中宋','fontsize',16);

subplot(1,3,3)
imshow(IF_image_phase)
set(gca,'fontname','times new roman','fontsize',16);    %设置图形对象属性
title('DOE相位分布','fontname','华文中宋','fontsize',16);

idea_image_bw=1 - im2bw(idea_image0, 0);
SSE = sum(sum((absF_image - idea_image).^2))/(sum(sum(idea_image.^2)))
n = sum(sum(((1-idea_image_bw).*absF_image).^2))/(sum(sum(absF_image.^2)))
n1 = sum(sum(((1-idea_image_bw).*absF_image).^2))/(sum(sum(idea_image.^2)))
n2 = sum(sum(((idea_image_bw).*absF_image).^2))/(sum(sum(idea_image.^2)))
sumnum = sum(sum(1-idea_image_bw));
mean_Abs_Fimage = sum(sum((1-idea_image_bw).*absF_image))/sumnum;
std = sqrt(sum(sum((1-idea_image_bw).*((absF_image - mean_Abs_Fimage).^2)))/(sumnum - 1))