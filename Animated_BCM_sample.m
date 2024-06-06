% Animated box-counting method calculation
% Using Kiran's data: Mitochondria fluorescence image
% 23 January 2024 - Devi Fitri Ferdania

clear
clc
close all
%% Load the data

% Counter 1 to 21

Arr = [];
for i = 1:21
    if i < 10
        counter = strcat("0",num2str(i));
    else
        counter = num2str(i);
    end
    Arr = [Arr counter];
end

%% Process of calculating BCM

% Calling the images
NameFile = imread("%076HCCCocl2Pax20024h.czi-HCCCocl2Pax20024h-0001.tif");

% Read the image
t = Tiff(NameFile,'r');
imageData = read(t);

% Specify some variables
m = size(imageData,1);
n = size(imageData,2);
s = max(m,n);
SizeSubMat = s/2;
epsilon = divisors(SizeSubMat);

% Create a binary image and add a block of black background
imbin = zeros(s,s);
for j=1:m
    for k=1:n
        if imageData(j,k) ~= 0
            imbin(j,k) = 1;
        end
    end
end

%%

% Count the boxes
Num = [];
Boxes = [];
Eps = [];
sz = size(epsilon,2);
figure
for j = sz:-1:1
    N = 0;
    Parts = s/epsilon(j);
    NumSubMat = (Parts)*(Parts);
    for k = 1:NumSubMat

        RowPos = ceil(k/Parts);
        ColPos = mod(k-1,Parts) + 1;

        r_st = ((RowPos - 1)*epsilon(j)) + 1;
        r_end = (RowPos*epsilon(j));

        c_st = ((ColPos - 1)*epsilon(j)) + 1;
        c_end = (ColPos*epsilon(j));

        % The submatrix
        SubMat = imbin(r_st:r_end,c_st:c_end);

        % Check whether a non zero matrix or not
        if all(SubMat(:)==0) == 0 % Not all entries are zero
            N = N+1;
        end
        
        % For animation
        h1 = subplot(1,2,1);
        imshow(imbin);
        axis on
        rectangle('Position',[c_st r_st epsilon(j) epsilon(j)],...
                'LineWidth',2,'LineStyle','--', 'EdgeColor', 'r')
        title(['Row = ', num2str(r_st),' Col = ', num2str(c_st), ' N = ', num2str(N)])
        hold on;
        
        if k == NumSubMat
            logN = log2(N);
            logEps = log2(1/epsilon(j));
            Num = [Num logN];
            Eps = [Eps logEps];
            plot(Eps, Num, 'r.')
            xlabel('log(1/epsilon)');
            ylabel('log(N)');
            xlim([-15 0]);
            ylim([0 30]);
            title(['Log-log plot for epsilon = ', num2str(epsilon(j))])
            hold on
        end


        if k == NumSubMat
            subplot(1,2,2);
            logN = log2(N);
            logEps = log2(1/epsilon(j));
            Num = [Num logN];
            Eps = [Eps logEps];
            plot(Eps, Num, 'r.')
            xlabel('log(1/epsilon)');
            ylabel('log(N)');
            xlim([-15 0]);
            ylim([0 30]);
            title(['Log-log plot for epsilon = ', num2str(epsilon(j))])
            hold on
        end

        MovieVector(k) = getframe;

        cla(h1); 
        xlabel('x')
        ylabel('y')

        pause(1);
    end
    Boxes = [Boxes N];
end

log_N = log2(Boxes);
l = size(log_N,2);
epsilons = zeros(1,l);
for j = l:-1:1
    epsilons(1,j) = 1/epsilon(1,(l+1)-j);
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
title(['Log-log plot with FD = ', num2str(FD)]);
hold on

