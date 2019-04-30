tic
close all;
clc;
clear all;
 
% Reading image as input
I = imread('lena.png');
Img = rgb2gray(I);
[m,n] = size(Img);
 
f=2;
%neighborhood window size = 2f+1 , ie 5x5
t=5;
% similarity window size = 2t+1 , ie 11x11
 
% Making gaussian kernel
su=1;          %standard deviation of gaussian kernel
sm=0;          % sum of all kernel elements (for normalization)
ks= 2*f+1;     % size of kernel (same as neighborhood window size)
% Initiating kernel
ker = zeros(ks,ks);
for x=1:ks
    for y=1:ks
        ab = x-f-1;   %horizontal distance of pixel from center(f+1, f+1)
        cd = y-f-1;   % vertical distance of pixel from center (f+1, f+1)
        ker(x,y) = 100*exp(((ab*ab)+(cd*cd))/(-2*(su*su)));
        sm = sm + ker(x,y);
    end
end
kernel = ker ./ f;
kernel = kernel / sm;   % normalization
 
% adding noise into image
noisex = imnoise(Img,'gaussian',0,0.002);
noisy = double(noisex);
 
% Assign a clear output image
cleared = zeros(m,n);
 
%Degree of filtering
h=10;
% Replicate boundaries of noisy image
noisy2 = padarray(noisy,[f,f],'symmetric');
 
% Now we'll calculate ouput for each pixel
for i=1:m
    for j=1:n
        im = i+f;   % to compensate for shift due to padarray function
        jn= j+f;
        % neighborhood of concerned pixel (we called it similarity window)
        W1 = noisy2(im-f:im+f , jn-f:jn+f);
        % BOundaries of similarity window for that pixel
        rmin = max(im-t, f+1);
        rmax = min(im+t, m+f);
        smin = max(jn-t, f+1);
        smax = min(jn+t, n+f);
        % Calculate weighted average next
        NL=0;    % same as cleared (i,j) but for simplicity
        Z =0;    % sum of all s(i,j)
        % Run loop through all the pixels in similarity window
        for r=rmin:rmax
            for s=smin:smax
                % neighborhood of pixel 'j' being compared for similarity
                W2 = noisy2(r-f:r+f, s-f:s+f);
                % square of weighted euclidian distances
                d2 = sum(sum(kernel.*(W1-W2).*(W1-W2)));
                % weight of similarity between both pixels : s(i,j)
                sij = exp(-d2/(h*h));
                % update Z and NL
                Z = Z + sij;
                NL = NL + (sij*noisy2(r,s));
            end
        end
        % normalization of NL
        cleared(i,j) = NL/Z;
    end
end
% convert cleared to uint8
cleared = uint8(cleared);
 
% show results
figure(1);
set(gcf, 'Position', get(0,'ScreenSize'));
subplot(1,2,1),imshow(noisex),title('noisy Image'),colormap(gray);
subplot(1,2,2),imshow(cleared),title('output of NL means filter'),colormap(gray);
toc