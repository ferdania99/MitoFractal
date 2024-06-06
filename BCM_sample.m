% This is a code to calculate fractal dimension
% of an image using the box counting method

clear all;
clc

%% Load the images
Im1 = imread("Kochsnowflake_402x464.jpg");

%% Pre-process the image
Im = Im1;
ImGray = rgb2gray(Im);
% imshow(ImGray)
ImBin = im2bw(ImGray, 0.9);
imshow(ImBin)

Convert 0 to 1 and 1 to 0
m = 500;
n = 400;

% maxSize = 1600;
ImMatrix = zeros(m,n);
for i = 1:464
    for j = 1:n
        if ImBin(i,j) == 0
            ImMatrix(i,j) = 1;
        end
    end
end

imshow(ImMatrix)

%% Start the box counting method
cdv = [1, 2, 4, 5, 10, 20, 25, 50, 100];

Num = [];
Boxes = [];
Eps = [];
sz = size(cdv,2);
for j = sz:-1:1
    N = 0;
    NumSubMat = (m/cdv(j))*(n/cdv(j));
    nbr = n/cdv(j);
    for k = 1:NumSubMat

        RowPos = ceil(k/nbr);
        ColPos = mod(k-1,nbr) + 1;

        r_st = ((RowPos - 1)*cdv(j)) + 1;
        r_end = (RowPos*cdv(j));

        c_st = ((ColPos - 1)*cdv(j)) + 1;
        c_end = (ColPos*cdv(j));

        % The submatrix
        SubMat = ImMatrix(r_st:r_end,c_st:c_end);

        % Check whether a non zero matrix or not
        if all(SubMat(:)==0) == 0 % Not all entries are zero
            N = N+1;
        end

        if k == NumSubMat
            % Uncomment the following line for animation
            % subplot(1,2,2);
            logN = log2(N);
            logEps = log2(1/cdv(j));
            Num = [Num logN];
            Eps = [Eps logEps];
            plot(Eps, Num, 'r.')
            xlabel('log(1/epsilon)');
            ylabel('log(N)');
            title('Log-log plot')
            hold on
        end
    end
    Boxes = [Boxes N];
end

Result2 = [fliplr(cdv)' Boxes'];

%% Box-counting dimension
log_N = log2(Boxes);
l = size(log_N,2);
epsilons = zeros(1,l);
for j = l:-1:1
    epsilons(1,j) = 1/cdv(1,(l+1)-j);
end
log_eps = log2(epsilons);

% Least Square

sigma_x = sum(log_eps);
sigma_y = sum(log_N);
sigma_xy = 0;
sigma_x2 = 0;
for j = 1:l
    sigma_xy = sigma_xy + (log_eps(1,j)*log_N(1,j));
    sigma_x2 = sigma_x2 + (log_eps(1,j))^2;
end

upper = (sz*sigma_xy) - (sigma_x*sigma_y);
lower = (sz*sigma_x2) - (sigma_x)^2;

FD = abs(upper/lower);

x = log_eps;
a = (sigma_y/l) - FD*(sigma_x/l);
y = FD*x + a;

figure
plot(log_eps,log_N,'ro');
hold on

plot(x,y,'b-');
title(['Log-log plot with box-counting dimension = ', num2str(FD)]);
xlabel('log(1/epsilon)');
ylabel('log(N)');
hold on
