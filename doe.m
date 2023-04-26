%% �������ͼ
image = im2double(rgb2gray(imread('matlab.jpg')));
image0 = imresize(image,[1024,1024],'bilinear');
idea_image0 = 1-image0;
figure(1)
subplot(1,3,1)
imshow(idea_image0);%��ʾ���������ǿ�ֲ�
set(gca,'fontname','times new roman','fontsize',16);    %����ͼ�ζ�������
title('��������۹�ǿ�ֲ�','fontname','��������','fontsize',16);


%% �����˹������������
Nxy = 1024;          %x,y���������
lambda = 632.8e-9;    %����Ϊ 632nm
k = 2*pi/lambda;    %����
w0 = 6e-3;        %��߳ߴ� 6mm
z = 0;
Z_R = pi*w0^2/lambda;           %��������
w_z = w0*sqrt(1+(z/Z_R)^2);     %������zλ�õİ뾶
%��������
[x,y] = meshgrid(linspace(-10e-3,10e-3,Nxy));     %����ֱ������
E = w0/w_z*exp(-(x.^2+y.^2)/w0^2)*exp(-1i*(k*z));%z=0����������ⳡ�ֲ�
I = (abs(E)).^2;%������������ǿ
phase = angle(E);%�������������λ
%% ����DOEԪ�������������λ
add_phase = (2*rand(Nxy,Nxy) - 1)*pi;
estimated_image = E.*exp(1i*add_phase); %ģ����ξ���DOE���˹����

%% ��������
tic
idea_image = sqrt(idea_image0);%�������ͼ
for m = 1:100
    F_image = fft2(estimated_image);%��k�ε���ʱ������ⳡ
    absF_image = abs(F_image)/max(max(abs(F_image)));
    sim = corrcoef(idea_image,absF_image);%��������ϵ��
    if sim>=0.98
        break
    end  
    F_image_phase = angle(F_image);%������λ
    F_estimated_image = idea_image.*exp(1i*F_image_phase);%������ⳡ�������������ⳡ������н���
    IF_image = ifft2(F_estimated_image);%������Ҷ�任�󣬼����k+1������ⳡ
    IF_image_phase = angle(IF_image);%������λ
    estimated_image = E.*exp(1i*IF_image_phase);%ʵ�ʵ�k+1������ⳡ
end
toc

subplot(1,3,2)
outimage = abs(F_image);
outimage = outimage/max(max(outimage));
imshow(outimage);
set(gca,'fontname','times new roman','fontsize',16);    %����ͼ�ζ�������
title('���������ֲ�','fontname','��������','fontsize',16);

subplot(1,3,3)
imshow(IF_image_phase)
set(gca,'fontname','times new roman','fontsize',16);    %����ͼ�ζ�������
title('DOE��λ�ֲ�','fontname','��������','fontsize',16);


idea_image6 = imresize(image,[128,128]);
idea_image6 = sqrt(idea_image6);
F_image6 = imresize(F_image,[128,128]);%kron(IF_image_phase - rand(64,64) * 0.01, ones(2));
F_image_phase6 = angle(F_image6);
F_estimated_image6 = idea_image6.*exp(1i*F_image_phase6);%������ⳡ�������������ⳡ������н���
IF_image6 = ifft2(F_estimated_image6);%������Ҷ�任�󣬼����k+1������ⳡ

IF_image_phase6 = angle(IF_image6);
[x,y] = meshgrid(linspace(-10e-3,10e-3,128));     %����ֱ������
E6 = w0/w_z*exp(-(x.^2+y.^2)/w0^2)*exp(-1i*(k*z));%z=0����������ⳡ�ֲ�
estimated_image6 = E6.*exp(1i*IF_image_phase6);%ʵ�ʵ�k+1������ⳡ
F_image6 = fft2(estimated_image6);
figure(2)
subplot(1,3,1)
imshow(idea_image6)

subplot(1,3,2)
outimage6 = abs(F_image6);
outimag6e = outimage6/max(max(outimage6));
imshow(outimag6e);%��ʾ���������ǿ�ֲ�

subplot(1,3,3)
imshow(IF_image_phase6)


