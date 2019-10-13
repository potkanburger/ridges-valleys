% parameters
clear all;
imgString = 'ellipses_generated';
imgExtension = '.jpg';
gray = true;

%program
fullImgPath = strcat(imgString,imgExtension);
ITEMP = imread(fullImgPath);
ITEMP = im2double(ITEMP);

global mathHp;

ITEMP = generateImg(ITEMP);
%ITEMP = double(rr);
%ITEMP = double((rr + (rr+15))/2);

if gray
    I = ITEMP;
else
    I = rgb2gray(ITEMP);
end
global Gx;
global Gy;
[Gx,Gy] = gaussfiltgradientxy(I,2);

global Hx;
global Hy;
[Hx, Hy] = hessianXgradient(I);

figure(1)
imshowpair(Gx,Gy,'montage')
title('Directional Gradients Gx and Gy, Using Central Difference')
saveas(gcf,strcat(imgString,'MATLAB.png'))

% pyXgrad = imread(strcat('x_gradient_',imgString,'.jpg'));
% pyYgrad = imread(strcat('y_gradient_',imgString,'.jpg'));
% 
% figure(1)
% imshowpair(pyXgrad,pyYgrad,'montage')
% title('Directional Gradients Gx and Gy, Using Central Difference, From Python')
% saveas(gcf,strcat(imgString,'PYT.png'))

mL = alternateheightPoints(I);

figure(2)
imshow(I);
hold on;
plotPoints(I, mL, 'r+');
title('HeightPoints');
saveas(gcf,strcat('hpts_',imgString,'MATLAB.png'))
hold off;

figure(3)
imshow(I);
hold on;
plotPointsBig(I, mL, 'ro');
title('HP with Gradient Vectors');
quiver(Gx, Gy);
quiver(Hx, Hy);

saveas(gcf,strcat('vectors_',imgString,'_MATLAB.png'))
hold off;

figure(4)
imshow(I);
hold on;

plotPointsBig(I, mL, 'ro');

clGx = onlyPx(mL, Gx);
clGy = onlyPx(mL, Gy);
clHx = onlyPx(mL, Hx);
clHy = onlyPx(mL, Hy);

title('Alternate HP with Gradient Vectors');
quiver(clGx, clGy);
quiver(clHx, clHy);
saveas(gcf,strcat('only_heightpoints_vectors_',imgString,'_MATLAB.png'))
hold off;

figure(5)
imshow(I);
hold on;
HPG1 = getHeightPoints(I, 1);
HPG15 = getHeightPoints(I, 1.5);
HPG2 = getHeightPoints(I, 2);
title('Height Points with Gaussians sigma 1, 1.5, 2');
plotPoints(I, HPG1, 'ro');
plotPoints(I, HPG15, 'bo');
plotPoints(I, HPG2, 'go');
hold off;

% [EVect1X, EVect1Y, EVect2X, EVect2Y] = eigvectorsHessian(I);
% 
% [HPG1, HP_NO_EIGVAL_CTRL] = getHeightPoints2(I, 1);
% 
% figure(6)
% imshow(I);
% hold on;
% quiver(Gx, Gy);
% quiver(EVect1X/2.0, EVect1Y/2.0, 'AutoScale','off');
% quiver(EVect2X/2.0, EVect2Y/2.0, 'AutoScale','off');
% plotPointsBig(I, HP_NO_EIGVAL_CTRL, 'go');
% plotPoints(I, HPG1, 'ro');
% hold off;


% figure(7)
% imshow(I);
% hold on;
% title('Height Points from Mathematical function');
% plotPoints(I, mathHp, 'r+');
% hold off;

% close all


function res = generateImg(inImg)
global mathHp;
%     res = zeros(size(inImg),'like',inImg);
%     [xSize, ySize] = size(inImg);
%     maxPoints = [100, 100];
%     maxPoints(:,:,2)= [230,230];
%     maxPoints = squeeze(maxPoints);
%     for x = 1:xSize
%         for y = 1:ySize
%             distance = xSize;
%             for p = 1:size(maxPoints,2)
%                 tmp_dst = sqrt((x - maxPoints(2,p))^2 + (y - maxPoints(1,p)/4)^2);
%                 if tmp_dst < distance
%                     distance = tmp_dst;
%                 end
%             end      
% 
%             val = double(255-distance);
%             res(y,x) = val;
%         end
%     end

alpha = pi*30/180;
a = 20;
b = 30;
sx = 0;
sy = 5;

sizeX = 50;
sizeY = 60;
ellipse = sqrt(((xx(sizeX,sizeY) - sx)*cos(alpha) + (yy(sizeX,sizeY) -sy)*sin(alpha))^2/a^2 + ((xx(sizeX,sizeY) - sx)*sin(alpha) - (yy(sizeX,sizeY) -sy)*cos(alpha))^2/b^2);

syms x;
syms y;
math_ellipse = sqrt(((x-sx)*cos(alpha) + (y-sy)*sin(alpha))^2/a^2 + ((x-sx)*sin(alpha) - (y-sy)*cos(alpha))^2/b^2);

gx = diff(math_ellipse, x);
gy = diff(math_ellipse, y);

gxx = diff(math_ellipse, x, 2);
gyy = diff(math_ellipse, y, 2);
gxy = diff(math_ellipse, x, y);

hessian_mat(1,1) = gxx;
hessian_mat(1,2) = gxy;
hessian_mat(2,1) = gxy;
hessian_mat(2,2) = gyy;

tt = [gx gy];
tt = tt.';
H = hessian_mat * tt;

hx = H(1);
hy = H(2);

determinant = hx*gy - hy*gx;
vals = 0;

for x = 1:sizeX
    previousSign = 0;
    previousDet = 0;
    for y = 1:sizeY
        try
            detL = double(subs(determinant));
            if y>1
                changeSign = sign(detL) ~= previousSign;
                
                if (almostZeroTolerance(detL, 0) || changeSign) && ~almostZero(double(subs(gx))) && ~almostZero(double(subs(gy))) 
                    
                    if(almostZeroTolerance(detL, 0))
                        valY = y;
                    else
                        if(abs(detL) <= abs(previousDet))
                            valY = y;
                        else
                            valY = y-1;
                        end
                    end
                    
                    mathHp(1,vals+1) = x;
                    mathHp(2,vals+1) = valY;
                    vals = vals+1;

                end                    
            end
            previousSign = sign(detL);
            previousDet = detL;
            
        catch
        end
    end
end

for y = 1:sizeY
    previousSign = 0;
    previousDet = 0;
    for x = 1:sizeX
        try
            detL = double(subs(determinant));
            if x>1
                changeSign = sign(detL) ~= previousSign;
                
                if (almostZeroTolerance(detL, 0) || changeSign) && ~almostZero(double(subs(gx))) && ~almostZero(double(subs(gy))) 
                    
                    if(almostZeroTolerance(detL, 0))
                        valX = x;
                    else
                        if(abs(detL) <= abs(previousDet))
                            valX = x;
                        else
                            valX = x-1;
                        end
                    end
                    
                    mathHp(1,vals+1) = valX;
                    mathHp(2,vals+1) = y;
                    vals = vals+1;

                end                    
            end
            previousSign = sign(detL);
            previousDet = detL;
            
        catch
        end
    end
end



res = double(ellipse);
end


function ptList = getHeightPoints(img, sigma)
    global Gx;
    global Gy;
    global Hx;
    global Hy;

    [Gx,Gy] = gaussfiltgradientxy(img,sigma);
    [Hx, Hy] = hessianXgradient(img);
    
    ptList = alternateheightPoints(img);
end

function [ptList, ptListSameEigVals] = getHeightPoints2(img, sigma)
    global Gx;
    global Gy;
    global Hx;
    global Hy;

    [Gx,Gy] = gaussfiltgradientxy(img,sigma);
    [Hx, Hy] = hessianXgradient(img);
    
    [ptList, ptListSameEigVals] = alternateheightPoints2(img);
end

function [determinant, tolerance] = getDeterminant(y, x)
    global Gx;
    global Gy;
    global Hx;
    global Hy;
    
    e = determinantTolerance();
    localH = [Hy(y,x) Hx(y,x)];
    
    htol = (1/norm(localH));    
    localH = localH/norm(localH);

    localHy = localH(1);
    localHx = localH(2);

    localG = [Gy(y,x) Gx(y,x)];
    
    gtol = (1/norm(localG)); 
    localG = localG/norm(localG);

    localGy = localG(1);
    localGx = localG(2);
    
    determinant = localHx*localGy - localHy*localGx;
    tolerance = abs((htol+e)*(gtol+e)) - abs((htol-e)*(gtol-e));
end



function res = determinantTolerance()
    res = 1e-38;
end

function res = almostZero(val)
    res = abs(val) <= 1e-38;   
end

function res = almostZeroTolerance(val, tolerance)
    res = abs(val) <= tolerance;   
end

function res = almostEqual(val1, val2)
    res = almostZeroTolerance(val1-val2, 1e-38) && sign(val1) == sign(val2);     
end

function [resX, resY] = gaussfiltgradientxy(I, sigma)
    It = imgaussfilt(I, sigma);
    [resX, resY] = imgradientxy(double(It), 'central');
end

function plotPoints(img, pts, ptMarker)

    if(size(pts)>0)
        xValues = pts(1,:);
        yValues = pts(2,:);
        plot(xValues, yValues, ptMarker, 'LineWidth', 2, 'MarkerSize', 2);
    end
end

function plotPointsBig(img, pts, ptMarker)
    if(size(pts)>0)
        xValues = pts(1,:);
        yValues = pts(2,:);
        plot(xValues, yValues, ptMarker , 'LineWidth', 4, 'MarkerSize', 4);
    end
end


function [gxx, gyy, gxy] = getSecondDerivatives(img)
  global Gx;
  global Gy;
  [gxx, gxy] = imgradientxy(double(Gx), 'central');
  [gyx, gyy] = imgradientxy(double(Gy), 'central');
end

function mat = hessian(gxx, gyy, gxy, x, y)
    mat(1,1) = gxx(y, x);
    mat(1,2) = gxy(y, x);
    mat(2,1) = gxy(y, x);
    mat(2,2) = gyy(y, x);
end

function [a,b,c] = getHessianABC(hessMat)
    a = hessMat(1,1);
    b = hessMat(1,2);
    c = hessMat(2,2);
end

function [eigvector1X, eigvector1Y, eigvector2X, eigvector2Y] = eigvectorsHessian(img)
    global Gx;
    global Gy;
    [ySize, xSize] = size(img);
    [iGxx, iGyy, iGxy] = getSecondDerivatives(img);
    eigvector1X = zeros(size(img),'like',img);
    eigvector1Y = zeros(size(img),'like',img);
    eigvector2X = zeros(size(img),'like',img);
    eigvector2Y = zeros(size(img),'like',img);
    
    for x = 1:xSize
        for y = 1:ySize
            localHessian = hessian(iGxx, iGyy, iGxy, x, y);
            [a,b,c] = getHessianABC(localHessian);
            [E,D] = eig(localHessian);
            eigvector1X(y, x) = E(1,1);
            eigvector1Y(y, x) = E(1,2);
            
            eigvector2X(y, x) = E(2,1);
            eigvector2Y(y, x) = E(2,2);
        end
    end

end
function [hessGradientX,hessGradientY] = hessianXgradient(img)
    global Gx;
    global Gy;
    [ySize, xSize] = size(img);
    [iGxx, iGyy, iGxy] = getSecondDerivatives(img);
    hessGradientX = zeros(size(img),'like',img);
    hessGradientY = zeros(size(img),'like',img);
    for x = 1:xSize
        for y = 1:ySize
            localHessian = hessian(iGxx, iGyy, iGxy, x, y);
            [a,b,c] = getHessianABC(localHessian);
            %cmp = [(a*Gx(y,x) + b*Gy(y,x)) (b*Gx(y,x) + c*Gy(y,x))];
            %hessGradientX(y, x) = a*Gx(y,x) + b*Gy(y,x);
            %hessGradientY(y, x) = b*Gx(y,x) + c*Gy(y,x);
            tt = [Gx(y,x) Gy(y,x)];
            tt = tt.';
            test = localHessian * tt;
            
            hessGradientX(y, x) = test(1);
            hessGradientY(y, x) = test(2);
        end
    end
end

function isLin = pixelIsLinear(hessVal, gxValue, gyValue)
    [a,b,c] = getHessianABC(hessVal);
    hessGradientX = a*gxValue + b*gyValue;
    hessGradientY = b*gxValue + c*gyValue;
    
    determinant = hessGradientX*gyValue - hessGradientY*gxValue;
    
    if (almostZero(hessGradientX) && almostZero(hessGradientY)) || (almostZero(gxValue) && almostZero(gyValue))
        isLin = false;
    else
        isLin = almostZeroTolerance(determinant, determinantTolerance());
    end
    %if isLin
    %    hessGradientX
    %    hessGradientY
    %    gxValue
    %   gyValue
    %end
    
end

function isLin = areLambdaLinear(hessMat, gxVal, gyVal)
    [a,b,c] = getHessianABC(hessMat);
    if ~almostZero(gyVal) && ~almostZero(gxVal)
        lambda1 = (a*gxVal + b*gyVal)/gxVal;
        lambda2 = (b*gxVal + c*gyVal)/gyVal;
        
        if almostEqual(lambda1, lambda2)
            isLin = true;
            return         
        end
        
    elseif (almostZero(gyVal) && ~almostZero(gxVal)) || (almostZero(gxVal) && ~almostZero(gyVal))
        isLin = almostZero(b*(gyVal+gxVal));
        return
    else
        isLin = true;
        return
    end
    isLin = false;
end


function cleared = onlyPx(ptList, mat)
    cleared = zeros(size(mat),'like',mat);
    [xSize, ySize] = size(ptList);
    [xSizeBis, ySizeBis] = size(mat);
    
    if ySize <= ySizeBis*xSizeBis
        for j = 1:ySize
            x = ptList(1,j);
            y = ptList(2,j);
            cleared(y,x) = mat(y,x);
        end
    end

end

function ptList = alternateheightPoints(img)
    global Gx;
    global Gy;
    
    global Hx;
    global Hy;
    
    [ySize, xSize] = size(img);
    [iGxx, iGyy, iGxy] = getSecondDerivatives(img);
    vals = 0;
    for x = 1:xSize
        previousSign = 0;
        previousDet = 0;
        for y = 1:ySize
            [determinant, tol] = getDeterminant(y,x);
            
            if y > 1
                changeSign = sign(determinant) ~= previousSign;
                if (almostZeroTolerance(determinant, tol) || changeSign) && (~almostZero(Hx(y,x)) || ~almostZero(Hy(y,x))) && (~almostZero(Gx(y,x)) || ~almostZero(Gy(y,x)))
                    hessVal = hessian(iGxx, iGyy, iGxy, x, y);
                    eigVals = eig(hessVal);
                    
                    if(almostZeroTolerance(determinant, tol))
                        valY = y;
                    else
                        if(abs(determinant) <= abs(previousDet))
                            valY = y;
                        else
                            valY = y-1;
                        end
                    end
                         
                    
                    if ~(almostZero(eigVals(1)) || almostZero(eigVals(2)) || almostEqual(eigVals(1), eigVals(2)))
                        ptList(1,vals+1) = x;
                        ptList(2,vals+1) = valY;
                        vals = vals+1;
                    end

                end
            end
            
            previousDet = determinant;
            previousSign = sign(determinant);
        end
    end
    
    for y = 1:ySize
        previousSign = 0;
        previousDet = 0;
        for x = 1:xSize
            [determinant, tol] = getDeterminant(y,x);
            
            if x > 1
                changeSign = sign(determinant) ~= previousSign;
                if (almostZeroTolerance(determinant, tol) || changeSign) && (~almostZero(Hx(y,x)) || ~almostZero(Hy(y,x))) && (~almostZero(Gx(y,x)) || ~almostZero(Gy(y,x)))
                    hessVal = hessian(iGxx, iGyy, iGxy, x, y);
                    eigVals = eig(hessVal);
                    
                    if(almostZeroTolerance(determinant, tol))
                        valX = x;
                    else
                        if(abs(determinant) <= abs(previousDet))
                            valX = x;
                        else
                            valX = x-1;
                        end
                    end
                         
                    
                    if ~(almostZero(eigVals(1)) || almostZero(eigVals(2)) || almostEqual(eigVals(1), eigVals(2)))
                        ptList(1,vals+1) = valX;
                        ptList(2,vals+1) = y;
                        vals = vals+1;
                    end

                end
            end
            
            previousDet = determinant;
            previousSign = sign(determinant);
        end
    end
    
    
    
    if(vals == 0)
        ptList = [];
    end
    
end


function [ptList, ptListSameEigVals] = alternateheightPoints2(img)
    global Gx;
    global Gy;
    
    global Hx;
    global Hy;
    fid=fopen("Log_heightpoints.txt", 'w');
    [ySize, xSize] = size(img);
    [iGxx, iGyy, iGxy] = getSecondDerivatives(img);
    vals = 0;
    valsE = 0;
    for x = 1:xSize
        for y = 1:ySize
            fprintf(fid, '[%i,%i]', x, y);
            
            [determinant, tol] = getDeterminant(y,x);
            fprintf(fid, 'determinant: %d, tolerance: %d, Hess: (%d, %d); Grad: (%d, %d)', determinant, tol, Hx(y,x), Hy(y,x), Gx(y,x), Gy(y,x));
            
            isHP = false;
            
            if almostZeroTolerance(determinant, tol) && (~almostZero(Hx(y,x)) || ~almostZero(Hy(y,x))) && (~almostZero(Gx(y,x)) || ~almostZero(Gy(y,x)))
                hessVal = hessian(iGxx, iGyy, iGxy, x, y);
                eigVals = eig(hessVal);
                
                if ~(almostZero(eigVals(1)) || almostZero(eigVals(2)))
                    ptListSameEigVals(1,vals+1) = x;
                    ptListSameEigVals(2,vals+1) = y;
                    valsE = valsE + 1;
                end
                
                if ~(almostZero(eigVals(1)) || almostZero(eigVals(2)) || almostEqual(eigVals(1), eigVals(2)))
                    ptList(1,vals+1) = x;
                    ptList(2,vals+1) = y;
                    vals = vals+1;
                    isHP = true;
                end
                fprintf(fid, ' ; eigenvalues : [%d,%d]', eigVals(1), eigVals(2));
                
            end
            
            if isHP
                fprintf(fid, ' ; is HeightPoint!\n');
            else
                fprintf(fid, ' ; is NOT heightpoint\n');
            end
        end
    end
    
    if(vals == 0)
        ptList = [];
    end
    
    if(valsE == 0)
        ptListSameEigVals = [];
    end
    
    fclose(fid);
end

function ptList = heightPoints(img)
    global Gx;
    global Gy;
    [ySize, xSize] = size(img);
    vals = 0;
    [iGxx, iGyy, iGxy] = getSecondDerivatives(img);
    for x = 1:xSize
        for y = 1:ySize
            gxValue = Gx(y,x);
            gyValue = Gy(y,x);       
            if ~almostZero(gxValue) || ~almostZero(gyValue)
                hessVal = hessian(iGxx, iGyy, iGxy, x, y);
                eigVals = eig(hessVal);
                if ~(almostZero(eigVals(1)) || almostZero(eigVals(2)) || almostEqual(eigVals(1), eigVals(2)))
                    if pixelIsLinear(hessVal, gxValue, gyValue)
                        ptList(1,vals+1) = x;
                        ptList(2,vals+1) = y;
                        vals = vals+1;
                    end
                end
            end
        end
    end
    
    if(vals == 0)
        ptList = [];
    end
end

