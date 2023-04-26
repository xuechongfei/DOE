%% �������ͼ
image = im2double(255 - rgb2gray(imread('1.jpg')));
image0 = imresize(image,[512,512],'bilinear');
idea_image0 = image0;
figure(2)
subplot(1,3,1)
imshow(idea_image0);%��ʾ���������ǿ�ֲ�
set(gca,'fontname','times new roman','fontsize',16);    %����ͼ�ζ�������
title('��������۹�ǿ�ֲ�','fontname','��������','fontsize',16);


%% �����˹������������
Nxy = 512;          %x,y���������
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
simlist = zeros(1,200);
idea_image_bw=1 - im2bw(idea_image0, 0);

P1 = idea_image/sum(sum(idea_image));
seta = 2;
for m = 1:200
    F_image = fft2(estimated_image);%��k�ε���ʱ������ⳡ
    absF_image = abs(F_image)*2/Nxy;
    sim = corrcoef(idea_image,absF_image);%��������ϵ��
    simlist(m) = sim(1,2);
    if sim(1,2)>=0.99
        break
    end  
    


    P2 = absF_image/sum(sum(absF_image));
    
    F_image_phase = angle(F_image);%������λ
    if m == 1
        F_image_phase1 = F_image_phase;
    end 
    
    leta = ((P1./P2).^seta);
    F_estimated_image = (idea_image + absF_image.*leta).*exp(1i*F_image_phase);%������ⳡ�������������ⳡ������н���
    IF_image = ifft2(F_estimated_image);%������Ҷ�任�󣬼����k+1������ⳡ
    IF_image_phase = angle(IF_image);%������λ
    estimated_image = E.*exp(1i*IF_image_phase);%ʵ�ʵ�k+1������ⳡ
end
toc

subplot(1,3,2)
outimage = abs(F_image);
outimage = outimage*2/Nxy;
% outimage = (outimage.^2);
imshow(outimage/max(max(outimage)));
set(gca,'fontname','times new roman','fontsize',16);    %����ͼ�ζ�������
title('���������ֲ�','fontname','��������','fontsize',16);

subplot(1,3,3)
imshow(IF_image_phase)
set(gca,'fontname','times new roman','fontsize',16);    %����ͼ�ζ�������
title('DOE��λ�ֲ�','fontname','��������','fontsize',16);

idea_image_bw=1 - im2bw(idea_image0, 0);
SSE = sum(sum((absF_image - idea_image).^2))/(sum(sum(idea_image.^2)))
n = sum(sum(((1-idea_image_bw).*absF_image).^2))/(sum(sum(absF_image.^2)))
n1 = sum(sum(((1-idea_image_bw).*absF_image).^2))/(sum(sum(idea_image.^2)))
n2 = sum(sum(((idea_image_bw).*absF_image).^2))/(sum(sum(idea_image.^2)))
sumnum = sum(sum(1-idea_image_bw));
mean_Abs_Fimage = sum(sum((1-idea_image_bw).*absF_image))/sumnum;
std = sqrt(sum(sum((1-idea_image_bw).*((absF_image - mean_Abs_Fimage).^2)))/(sumnum - 1))