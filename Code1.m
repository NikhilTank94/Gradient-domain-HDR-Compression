%%Matlab code for Gradient Domain HDR compression
% Sequence
% taking log of luminance
% taking gradient of log(luminance)
% obtaning the phi
% cal gradient attenuation function for each level of gauss
% scaling factor at every level
% cal gradient attenuation function for each level of gauss
% obtaning cap Phi
% multiplying gradx,y to phi and solving poisson equation
% equating lap I = div G
% make verticle source matrix
% convert verticle into required
% transform the HDR color according to the new calculated luminance
clc
clear all
count=0;
hdr = hdrread('H.hdr');
hdr = imresize(hdr, [2*64 2*128]);
lab = rgb2lab(hdr);
L_k = lab(:, :, 1); % Extract the L image
A = lab(:, :, 2);   % Extract the a image
B = lab(:, :, 3);   % Extract the b image
GDHDR = hdr;
% taking log of luminance
lnL=abs(log(L_k));

% taking gradient of log(luminance)
[gradx,grady] = gradient(lnL);

%% obtaning the phi
% make gauss pyr till d
n = log2(size(L_k));
d = min(n)-5;   % to make course matrix with min height or width = 32 
G = lnL;
gauss = {G};
for i = 1:d 
    G = impyramid(G,'reduce');
    gauss = [gauss;G];
end

% cal gradient attenuation function for each level of gauss
gradH = {};
for k = 1:d+1
    temp = gauss{k};
    [x,y] = size(temp);
    value = zeros(x, y, 2, 'double');
        for i = 1:x
            for j = 1:y
                if i+1>x
                    value(i,j,1) = (0-temp(i-1,j))/(2^(k+1));
                    %value(i,j,1) = (0)/(2^(k+1));
                elseif i-1<1
                    value(i,j,1) = (temp(i+1,j)-0)/(2^(k+1));
                    %value(i,j,1) = (0)/(2^(k+1));
                else    
                    value(i,j,1) = (temp(i+1,j)-temp(i-1,j))/(2^(k+1)); 
                end
                if j+1>y
                    value(i,j,2) = (0-temp(i,j-1))/(2^(k+1));
                    %value(i,j,2) = (0)/(2^(k+1));
                elseif j-1<1
                    value(i,j,2) = (temp(i,j+1)-0)/(2^(k+1));
                    %value(i,j,2) = (0)/(2^(k+1));
                else    
                    value(i,j,2) = (temp(i,j+1)-temp(i,j-1))/(2^(k+1));
                end
            end
        end
        gradH = [gradH;value];
end

% scaling factor at every level
% cal gradient attenuation function for each level of gauss
phi = {};
for k = 1:d+1
    temp1 = gradH{k};
    [x,y,z] = size(temp1);
    value = zeros(x, y, 'double');
    alpha = 0.1*(mean2(abs(temp1(:,:,1))) + mean2(abs(temp1(:,:,2))) )/2;    
    beta = 0.88;
    for i = 1:x
            for j = 1:y
                norm_gradH = sqrt(temp1(i,j,1)^2 + temp1(i,j,2)^2);
                value(i,j) = (alpha/norm_gradH)*(norm_gradH/alpha).^beta;
            end
    end
    phi = [phi;value];
end

% obtaning cap Phi 
 P=phi{d+1};
for k = 1:d
     app = phi{d+1-k};
     org = imresize(P,2);
     [x,y] = size(app);
     for i = 1:x
         for j = 1:y
             P(i,j) = org(i,j)*app(i,j);
         end
     end
end

%% multiplying gradx,y to phi and solving poisson equation
[x,y] = size(P);     
Gx = zeros(x, y, 'double');
for i = 1:x
    for j = 1:y
             Gx(i,j) = gradx(i,j)*P(i,j);
    end
end

Gy = zeros(x, y, 'double');
for i = 1:x
    for j = 1:y
             Gy(i,j) = grady(i,j)*P(i,j);
    end
end



%% equating lap I = div G 
[x,y] = size(Gx);
divG = zeros(x, y, 'double');
for i = 1:x
    for j = 1:y
        if i-1<1
            divG(i,j) = 0;
        elseif j-1<1
            divG(i,j) = 0; 
        else
            divG(i,j) = Gx(i,j) - Gx(i-1,j) + Gy(i,j) - Gy(i,j-1);
        end
    end     
end

poi = speye(x*y,x*y);
for i =1:x*y %row variable %start & endpt included
    poi(i,i) = (-4);
    if i-y >= 1
        poi(i,i-y) = 1;
    end
    if i-1 >= 1
        if rem((i-1),y) ~= 0
            poi(i,i-1) = 1;
        end
    end
    if i+y<=x*y
        poi(i,i+y) = 1;
    end
    if i+1<=x*y
        if rem(i,y) ~= 0
            poi(i,i+1) = 1;
        end
    end
end

%% make verticle source matrix
[x,y] = size(Gx);
b = zeros(x*y,1);
r=0;
for i = 1:x %row variable
    for j = 1:y
        b(i+j+r-1,1) = divG (i,j);
    end
    r = r + y - 1;
end

vertI=poi\b;

%% convert verticle into required
I = zeros(x,y);
r=0;
for i = 1:x %row variable
    for j = 1:y
        I(i,j) = vertI(i+j+r-1,1) ;
    end
    r = r + y - 1;
end
 L_j = exp(I);
lab(:, :, 1) = exp(I);
HDR = lab2rgb(lab);
[x,y,z] = size(hdr);
s=.25;
GDHDRC = zeros(x,y,z,'double');

%% transform the HDR color according to the new calculated luminance
for k=1:z
    for i = 1:x
         for j = 1:y                     
             GDHDR(i,j,k) = L_j(i,j).*(hdr(i,j,k)/L_k(i,j)).^s;                          
         end
    end
end
imwrite(GDHDR,'Compressed_HDR.png'); %Gradient domain compressed HDR image
imwrite(HDR,'HDR.png'); %untoned hdr image
