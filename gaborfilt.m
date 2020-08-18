%%Mohammadreza Pourfard (PhD student)
% 20 day 1392 
% Amirkabir university of Iran

% 1-This file get an image and convolve it with a bank of Gabor filter 
% with different orientations
% 2- This file compute the maximum variance of gabor filter result for
% every block and construct a fuuzzy color maximum variance gabor filter image
% 3- This file shows the angular histogram of gabor filter result and smooth
% the histogram and cut it below a threshold and find its peaks
% 4-This file unsupervisedly segment the colored maximum variance gabor
% filter image and shows and saves the final result

clear all;
close all;
clc;
% I=imread('C:\Users\pourfardm\Pictures\barbara (1).png');
% I=rgb2gray(imread('C:\Users\pourfardm\Pictures\Shingu01.bmp'));%No rotated region
% I=rgb2gray(imread('C:\Users\pourfardm\Pictures\Shingu02-01-45.bmp'));%one rotated region
% I=rgb2gray(imread('C:\Users\pourfardm\Pictures\Shingu02-04-315.tif'));%three rotated region
I=rgb2gray(imread('C:\Users\HP\Desktop\Extract\nano code- Gabor filter\original pics.jpg'));%special pics.jpg

% I=rgb2gray(imread('C:\Users\office\Documents\MATLAB\angles_nano\1.tif'));%No rotated region
% I=rgb2gray(imread('C:\Users\office\Documents\MATLAB\angles_nano\2.tif'));%one rotated region
% I=rgb2gray(imread('D:\MATLAB_D\angles_nano\3.tif'));%three rotated region
% I=rgb2gray(imread('D:\MATLAB_D\angles_nano\angles_nano\4.tif'));%special pics.jpg
I=imresize(I,0.5,'box');

% I=imrotate(I,90);

% [I2,bw,bw2,label,x,y]=adaptivethresh3(I);
bw2=I/255;
% imshow(I);
J=double(bw2);
% J=imresize(J,0.5);

%% Crop
% I2 = imcrop(I);
% figure, imshow(I2)
theta = 3;
sigma = 0.65*theta;
%% main parameter
filterSize = 20;
Range_size=8;
% last_thetha=0.34;
last_thetha=91;
[v1 v2]=meshgrid(1:filterSize,1:filterSize);
%%
iteration=0;
for zav=1:1:last_thetha
    iteration=iteration+1;%Proportional to angle
    phi =zav/180*pi;
    phi_all(iteration)=zav;
%     G = zeros(filterSize);%Gabor filter
    
    for i=(0:filterSize-1)/filterSize
        for j=(0:filterSize-1)/filterSize
            xprime= j*cos(phi);
            yprime= i*sin(phi);
            K = exp(2*pi*theta*sqrt(-1)*(xprime+ yprime));
            G(round((i)*filterSize)+1,round((j)*filterSize)+1) = exp(-(i^2+j^2)/(sigma^2))*K;
        end
        
    end
    
    K = conv2(J,G);
    U(:,:,iteration)=imag(K);
%     colormap(gray);
figure(1);
    subplot(2,2,1),imshow(I);title('Original image');
    subplot(2,2,2),imshow(J);title('Scaled (0-255) Original image');
    subplot(2,2,3),imagesc(U(:,:,iteration));title(['Gabor filter bank results: Angle is: ',num2str(zav)]);axis square;
    %     subplot(1,2,2),imagesc(U>0*max(U(:)));title(['Angle is: ',num2str(angle*180/8)]);
%     imshow(U(:,:,iteration));
    subplot(2,2,4),imshow(imag(G));title('Gabor filter');
%     colormap(jet);%color
    figure(2),mesh(v1,v2,imag(G));title(['Mesh: Angle is: ',num2str(zav)]);axis square;
    xlabel('X');ylabel('Y');zlabel('Gabor');
%     CC=imagesc(U(:,:,iteration));
%     subplot(2,3,6),imshow(imread(CC));
    
%     pause(1);
    drawnow;
    HH=imag(G);
end
% [Index value] = max(U,[],3);
% figure,imshow(value);
%%


sz=size(I);
q1=ceil(sz(1)/Range_size)*Range_size;
q2=ceil(sz(2)/Range_size)*Range_size;

for i=1:iteration %size(U,3)=number of angles
    UU=U(:,:,i);
    II=padarray(UU,[q1-sz(1),q2-sz(2)],'post');
    
    SZ=size(II);
    %put image between 0-255
%     xmin=min(II(:));
%     xmax=max(II(:));
%     a=0;
%     b=255;
%     III=(b-a)*(II(:)-xmin)/(xmax-xmin)+a;
%     II=reshape(III,[SZ(1) ,SZ(2)]);
    
    II=double(II)/255; %Normalize II
    R(:,:,i)=im2col(II,[Range_size,Range_size],'distinct');  %blocks of image
    [i iteration]
end

E=var(R,1);
EE=shiftdim(E);
% EE = col2im(E,[1 1],[SZ(1) SZ(2)],'distinct');
% EE=reshape(E,ceil([SZ(1)/Range_size SZ(2)/Range_size size(R,3)]));
% window=[0.1 0.4 0.7 1 0.7 0.4 0.1 ];
window=[0.7 0.8 0.9 1 0.9 0.8 0.7];
window=window/sum(window(:));
%% Smoothing the variance signal

for i=1:size(EE,1)
    de=conv(EE(i,:),window);
%     de=corr(EE(i,:),window);
    FF(i,:)=de(floor(size(window,2)/2):size(de,2)-floor(size(window,2)/2)-1)';
end

% [bbx bby]=meshgrid(1:size(EE,2),1:size(EE,1));
% surf(bbx,bby,EE);title('Surf original');




% save gaborimage.mat I J U R phi_all iteration;

% BB=shiftdim(mean(R,1));
% 


%% reconstruct
%find the maximum correlation block between all angles
% [RR index]=max(var(R,1),[],3); %R=blocks of image
%or
[RR index]=max(FF,[],2); %R=blocks of image
RR=RR';index=index';

sz=size(II);
m=1;
% zaviyeh=90;

Kr=m*trimf(0:iteration-1,[0,.25,1]*iteration);
Kg=m*trimf(0:iteration-1,[0,.75,1]*iteration);
Kb=m*max(trimf(0:iteration-1,[-1,0,.5]*iteration),...
    trimf(0:iteration-1,[.5,1,2]*iteration));% for two part trimf

for j=1:size(index,2)%number of blocks
    Q(:,j)=R(:,j,index(j));%Paste the best matched block
    color(1:Range_size^2,j)=index(j);%index of the best match
end
I_final=col2im(Q,[Range_size Range_size],[sz(1) sz(2)],'distinct');
I_color=col2im(color,[Range_size Range_size],[sz(1) sz(2)],'distinct');

gamma=0.5;
I_final2(:,:,1)=gamma*I_final+(1-gamma)*Kr(I_color);
I_final2(:,:,1)=I_final2(:,:,1)/max(max(I_final2(:,:,1)));
I_final2(:,:,2)=gamma*I_final+(1-gamma)*Kg(I_color);
I_final2(:,:,2)=I_final2(:,:,2)/max(max(I_final2(:,:,2)));
I_final2(:,:,3)=gamma*I_final+(1-gamma)*Kb(I_color);
I_final2(:,:,3)=I_final2(:,:,3)/max(max(I_final2(:,:,3)));

% figure;
% subplot(2,2,1),imshow(uint8(I_final*255));title('Selected region');
% subplot(2,2,2),imshow(I_final2);title(['Color Selected region with block size: ',num2str(Range_size)]);
% subplot(2,2,3),imshow(J);title('Scaled (0-255) Original image');
%% histogram of angles
angles=zeros(1,iteration);
for i=1:size(index,2)
    angles(index(i))=angles(index(i))+1;
end
figure,plot(angles);title('Histogram of angles of final constructed image');xlabel('Angle with degree');
hold on;

n2=10;
filt2=hamming(n2);filt2=filt2/sum(filt2);
angles1=conv(angles,filt2);
angles1=circshift(angles1,[0 round(-n2/2)+1]);
angles1=angles1(1:iteration);
plot(angles1,'r');

angles2=angles;
threshold=0.2;
angles2(angles<max(angles)*threshold)=0;
plot(angles2,'g');legend('Normal histogram','Smoothed histogram','Thresholded Smoothed histogram');
%%%%%%%%%%%%%%%%%%%%%%%
%% unsupervised segmentation
threshold=0.2;

%the k largest value of angles1
% k=5;
% sorted_angles1=sort(angles1,'descend');%sorted angles1
% angles4=zeros(size(angles1));
% angles4(angles1>=sorted_angles1(k))=angles1(angles1>=sorted_angles1(k));

H2(1:last_thetha)=angles;%value of each angle**************************
L1=H2;
L2=circshift(H2,[0 1]);
L3=circshift(H2,[0 -1]);
L4=circshift(H2,[0 2]);
L5=circshift(H2,[0 -2]);
% L6=circshift(H2,[0 3]);
% L7=circshift(H2,[0 -3]);
ind=find((L1>L2) & (L1>L3) & (L1>L4) & (L1>L5));% &(L1>L6) & (L1>L7));
[sorted_L IX]=sort(L1(ind),'descend');%sorted peaks

k_threshold=1;
% k=3;%**********************************************************************
k=round(numel(ind)*k_threshold); %the k largest value of angles1

num_peaks=sum((L1(ind)>threshold*sorted_L(1)));
plot(ind,L1(ind),'c*');
plot(ind,(L1(ind)>threshold*sorted_L(1)) .*L1(ind),'m.');


for iter=1:iteration
%     [v INDEX]=min((iter-ind).^2);
    [v INDEX]=min((iter-ind(IX(1:k))).^2);
%     unsupervised_region(iter)=L1(ind(INDEX));
%     unsupervised_region_label(iter)=ind(INDEX);
    unsupervised_region(iter)=L1(ind(IX(INDEX)));
    unsupervised_region_label(iter)=ind(IX(INDEX));
end
plot(unsupervised_region,'y');
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% iteration1=numel(ind);
% Kr=m*trimf(0:iteration1-1,[0,.25,1]*iteration1);
% Kg=m*trimf(0:iteration1-1,[0,.75,1]*iteration1);
% Kb=m*max(trimf(0:iteration1-1,[-1,0,.5]*iteration1),...
%     trimf(0:iteration1-1,[.5,1,2]*iteration1));% for two part trimf



gamma=0.5;
I_final3(:,:,1)=gamma*I_final+(1-gamma)*Kr(unsupervised_region_label(I_color));
I_final3(:,:,1)=I_final3(:,:,1)/max(max(I_final3(:,:,1)));
I_final3(:,:,2)=gamma*I_final+(1-gamma)*Kg(unsupervised_region_label(I_color));
I_final3(:,:,2)=I_final3(:,:,2)/max(max(I_final3(:,:,2)));
I_final3(:,:,3)=gamma*I_final+(1-gamma)*Kb(unsupervised_region_label(I_color));
I_final3(:,:,3)=I_final3(:,:,3)/max(max(I_final3(:,:,3)));

I_color2=unsupervised_region_label(I_color);

%% color real image
m=1;
xt2=0;%number of domains after unsupervised segmentation
for i=1:iteration
    x2=sum(find(I_color2==i));
    if(x2>0)
        xt2=xt2+1;
        x2_index(xt2)=i;
    end
end
% zaviyeh=90;
% iteration=xt2;
% Krr=m*trimf(0:iteration-1,[0,.25,1]*iteration);
% Kgg=m*trimf(0:iteration-1,[0,.75,1]*iteration);
% Kbb=m*max(trimf(0:iteration-1,[-1,0,.5]*iteration),...
%     trimf(0:iteration-1,[.5,1,2]*iteration));% for two part trimf
I_final4(:,:,1)=double(I).*Kr(I_color2(1:size(I,1),1:size(I,2)));
I_final4(:,:,2)=double(I).*Kg(I_color2(1:size(I,1),1:size(I,2)));
I_final4(:,:,3)=double(I).*Kb(I_color2(1:size(I,1),1:size(I,2)));



figure;%,subplot(2,3,1)
subplot(2,3,2),imshow(uint8(I_final*255));title('Selected region');
subplot(2,3,3),imshow(I_final2);title(['Color Selected region with block size: ',num2str(Range_size),' with: ' ,num2str(numel(L1~=0)),' Color Regions']);
subplot(2,3,4),imshow(J);title('Scaled (0-255) Original image');
% subplot(2,3,5),imshow(I_final3);title(['Color unsupervised segmentation with: ',num2str(numel(ind)),' Regions']);
subplot(2,3,5),imshow(I_final3);title(['Color unsupervised segmentation with: ',num2str(k),' Color Regions']);
subplot(2,3,6),imshow(uint8(I_final4));title('Color unsupervised segmentation in real image');

save gaborimage.mat I II J U R phi_all iteration Kr Kg Kb I_final I_final2 I_final3 I_color last_thetha angles angles1 angles2 last_thetha k unsupervised_region_label unsupervised_region ind IX Range_size last_thetha threshold k_threshold I_color2;
save gaborcolor.mat I_color I_color2 I Range_size;
save Kbig.mat Kr Kg Kb;





