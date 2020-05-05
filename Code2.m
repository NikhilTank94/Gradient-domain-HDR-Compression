%%Matlab code for Image matting
clc
clear all
count=0;
gre = imread('Green.png'); %RGB order of color
gre = imresize(gre, .5);
blu = imread('Blue.png');
blu = imresize(blu, .5);
BG_new = imread('BG.png');
BG_new = imresize(BG_new, .5); % new backgroung
[x,y,z]=size(gre);
alpha = zeros(x, y);
FG = gre;
IMG = BG_new;
a = (blu-gre);
b = (gre-blu);
c = 254;  %c = (gre(1,1,2)-blu(1,1,2))
for k = 1:z % calc alpha only for green channel
    for i = 1:x
        for j = 1:y 
            if k ==2
                alpha(i,j,k) = (254-b(i,j,k));
            elseif k ==3
                alpha(i,j,k) = (254-a(i,j,k));
            end
        end
    end
end
Alpha = alpha./c;
Alpha(:,:,1) = Alpha (:,:,2);
Beta = ones(x,y,z) - Alpha; %(1 - Alpha(i,j,k)) at every pt
BG = zeros(x, y, z);
BG(:,:,2) = 254;
for k = 1:z % calc alpha only for green channel
    for i = 1:x
        for j = 1:y  
            FG(i,j,k) = (gre(i,j,k)-Beta(i,j,k).*BG(i,j,k))./Alpha(i,j,k);
        end
    end
end
for k = 1:z % ransfer the same effect to a new back ground of choice
    for i = 1:x
        for j = 1:y  
            IMG(i,j,k) = Alpha(i,j,k).*FG(i,j,k) + Beta(i,j,k).*BG_new(i,j,k);
        end
    end
end
imwrite(FG,'FG.png');
imwrite(IMG,'Image_matted.png');

