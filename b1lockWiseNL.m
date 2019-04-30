tic
close all;
clc;
clear all;
 
% Reading image as input
I = imread('lena.png');
Img = rgb2gray(I);
[m,n] = size(Img);
 
f=2;
%neighborhood window size = 2f+1 
t=7;
% similarity window size = 2t+1 
 
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
noisy2 = padarray(noisy,[4,4],'symmetric');
 
% Now we'll calculate ouput for each pixel
for i=1:m
    for j=1:n
        im = i+4;   % to compensate for shift due to padarray function
        jn= j+4;
        % neighborhood of concerned pixel (we called it similarity window)

        % initialize 4 different reference pixels.      
        p1x = im - f;
        p1y = jn - f;
        p2x = im + f;
        p2y = jn - f;
        p3x = im - f;
        p3y = jn + f;
        p4x = im + f;
        p4y = jn + f;
        wa = noisy2(p1x-2:p1x+2 , p1y-2:p1y+2);
        wb = noisy2(p2x-2:p2x+2 , p2y-2:p2y+2);
        wc = noisy2(p3x-2:p3x+2 , p3y-2:p3y+2);
        wd = noisy2(p4x-2:p4x+2 , p4y-2:p4y+2);
        % initialize the volumn for the corresponding blocks
        p1rmin = max(p1x-t, 5);
        p1rmax = min(p1x+t, m+4);
        p1smin = max(p1y-t, 5);
        p1smax = min(p1y+t, n+4);
        
        p2rmin = max(p2x-t, 5);
        p2rmax = min(p2x+t, m+4);
        p2smin = max(p2y-t, 5);
        p2smax = min(p2y+t, n+4);
        
        p3rmin = max(p3x-t, 5);
        p3rmax = min(p3x+t, m+4);
        p3smin = max(p3y-t, 5);
        p3smax = min(p3y+t, n+4);
        
        p4rmin = max(p4x-t, 5);
        p4rmax = min(p4x+t, m+4);
        p4smin = max(p4y-t, 5);
        p4smax = min(p4y+t, n+4);
        
        % Calculate weighted average
        NL1=0;    % same as cleared (i,j) but for simplicity
        Z1 =0;    % sum of all w(i,j)
        % Run loop through all the pixels in similarity window
        for r=(p1rmin+2):5:p1rmax
            for s=(p1smin+2):5:p1smax
                % neighborhood of pixel 'j' being compared for similarity
                W2 = noisy2(r-f:r+f, s-f:s+f);
                % square of weighted euclidian distances
                d2 = sum(sum(kernel.*(wa-W2).*(wa-W2)));
                % weight of similarity between both pixels : s(i,j)
                wij = exp(-d2/(h*h));
                % update Z and NL
                Z1 = Z1+ wij;
                NL1 = NL1 + (wij*noisy2(r,s));
            end
        end
        % normalization of NL
        NL1 = NL1/Z1;
      
                % Calculate weighted average next
        NL2=0;    % same as cleared (i,j) but for simplicity
        Z2 =0;    % sum of all w(i,j)
        % Run loop through all the pixels in similarity window
        for r=(p2rmin+2):5:p2rmax
            for s=(p2smin+2):5:p2smax
                % neighborhood of pixel 'j' being compared for similarity
                W2 = noisy2(r-f:r+f, s-f:s+f);
                % square of weighted euclidian distances
                d2 = sum(sum(kernel.*(wb-W2).*(wb-W2)));
                % weight of similarity between both pixels : s(i,j)
                wij = exp(-d2/(h*h));
                % update Z and NL
                Z2 = Z2+ wij;
                NL2 = NL2 + (wij*noisy2(r,s));
            end
        end
        % normalization of NL
        NL2 = NL2/Z2;
        
        
                      
        NL3=0;    % same as cleared (i,j) but for simplicity
        Z3 =0;    % sum of all w(i,j)
        % Run loop through all the pixels in similarity window
        for r=(p3rmin+2):5:p3rmax
            for s=(p3smin+2):5:p3smax               
                W2 = noisy2(r-f:r+f, s-f:s+f);
                % square of weighted euclidian distances
                d2 = sum(sum(kernel.*(wc-W2).*(wc-W2)));
                % weight of similarity between both pixels : s(i,j)
                wij = exp(-d2/(h*h));
                % update Z and NL
                Z3 = Z3+ wij;
                NL3 = NL3 + (wij*noisy2(r,s));
            end
        end
        % normalization of NL
        NL3 = NL3/Z3;
        
        
                      
        NL4=0;    % same as cleared (i,j) but for simplicity
        Z4 =0;    % sum of all w(i,j)
        % Run loop through all the pixels in similarity window
        for r=(p4rmin+2):5:p4rmax
            for s=(p4smin+2):5:p4smax
                % neighborhood of pixel 'j' being compared for similarity
                W2 = noisy2(r-f:r+f, s-f:s+f);
                % square of weighted euclidian distances
                d2 = sum(sum(kernel.*(wd-W2).*(wd-W2)));
                % weight of similarity between both pixels : s(i,j)
                wij = exp(-d2/(h*h));
                % update Z and NL
                Z4 = Z4+ wij;
                NL4 = NL4 + (wij*noisy2(r,s));
            end
        end
        % normalization of NL
        NL4 = NL4/Z4;
        % new images generated
        cleared(i,j) = (NL1 + NL2 + NL3 + NL4)/4;
    
    end
end
% convert cleared to uint8
cleared = uint8(cleared);
 
% show results
figure(1);
set(gcf, 'Position', get(0,'ScreenSize'));
subplot(1,2,1),imshow(noisex),title('noisy Image'),colormap(gray);
subplot(1,2,2),imshow(cleared),title('output of blockwise NL means filter'),colormap(gray);
toc