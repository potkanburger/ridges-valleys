% parameters

imgString = 'ellipses_generated';
imgExtension = '.jpg';
gray = true;

%program
fullImgPath = strcat(imgString,imgExtension);
ITEMP = imread(fullImgPath);
ITEMP = im2double(ITEMP);
if gray
    I = ITEMP;
else
    I = rgb2gray(ITEMP);
end
global Gx;
global Gy;
[Gx,Gy] = gaussfiltgradientxy(I,2);

figure(1)
imshowpair(Gx,Gy,'montage')
title('Directional Gradients Gx and Gy, Using Central Difference')
saveas(gcf,strcat(imgString,'MATLAB.png'))

pyXgrad = imread(strcat('x_gradient_',imgString,'.jpg'));
pyYgrad = imread(strcat('y_gradient_',imgString,'.jpg'));
figure(2)

imshowpair(pyXgrad,pyYgrad,'montage')
title('Directional Gradients Gx and Gy, Using Central Difference, From Python')
saveas(gcf,strcat(imgString,'PYT.png'))

mL = heightPoints(I);

figure(3)
imshow(I);
hold on;
plotPoints(I, mL, 'r+');
title('HeightPoints');
saveas(gcf,strcat('hpts_',imgString,'MATLAB.png'))
hold off;

global Hx;
global Hy;
[Hx, Hy] = hessianXgradient(I);

figure(4)
imshow(I);
hold on;
plotPointsBig(I, mL, 'ro');
title('HP with Gradient Vectors');
quiver(Gx, Gy);
quiver(Hx, Hy);

saveas(gcf,strcat('vectors_',imgString,'_MATLAB.png'))
hold off;

figure(5)
imshow(I);
hold on;
amL = alternateheightPoints(I);

plotPointsBig(I, amL, 'ro');

clGx = onlyPx(amL, Gx);
clGy = onlyPx(amL, Gy);
clHx = onlyPx(amL, Hx);
clHy = onlyPx(amL, Hy);

title('Alternate HP with Gradient Vectors');
quiver(clGx, clGy);
quiver(clHx, clHy);
saveas(gcf,strcat('only_heightpoints_vectors_',imgString,'_MATLAB.png'))
hold off;

figure(6)
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

[EVect1X, EVect1Y, EVect2X, EVect2Y] = eigvectorsHessian(I);

figure(7)
imshow(I);
hold on;
quiver(Gx, Gy);
quiver(EVect1X/10.0, EVect1Y/10.0, 'AutoScale','off');
quiver(EVect2X/10.0, EVect2Y/10.0, 'AutoScale','off');
plotPoints(I, HPG1, 'ro');
hold off;
% close all

function ptList = getHeightPoints(img, sigma)
    global Gx;
    global Gy;
    global Hx;
    global Hy;

    [Gx,Gy] = gaussfiltgradientxy(img,sigma);
    [Hx, Hy] = hessianXgradient(img);
    
    ptList = alternateheightPoints(img);
end

function res = almostZero(val)
    res = abs(val) <= 1e-38;   
end

function res = almostEqual(val1, val2)
    res = almostZero(val1-val2);     
end

function [resX, resY] = gaussfiltgradientxy2(I, sigma)
    [tmpX, tmpY] = imgradientxy(double(I), 'central');
    resX = imgaussfilt(tmpX, sigma);
    resY = imgaussfilt(tmpY, sigma);
end

function [resX, resY] = gaussfiltgradientxy(I, sigma)
    It = imgaussfilt(I, sigma);
    [resX, resY] = imgradientxy(double(It), 'central');
end

function plotPoints(img, pts, ptMarker)
    xValues = pts(1,:);
    yValues = pts(2,:);
    plot(xValues, yValues, ptMarker, 'LineWidth', 2, 'MarkerSize', 2);
end

function plotPointsBig(img, pts, ptMarker)
    xValues = pts(1,:);
    yValues = pts(2,:);
    plot(xValues, yValues, ptMarker , 'LineWidth', 4, 'MarkerSize', 4);
end


function [gxx, gyy, gxy] = getSecondDerivatives(img)
  global Gx;
  global Gy;
  [gxx, gxy] = imgradientxy(double(Gx), 'central');
  [gyx, gyy] = imgradientxy(double(Gy), 'central');
end

function mat = hessian(gxx, gyy, gxy, x, y)
    mat(1,1) = gxx(x,y);
    mat(1,2) = gxy(x,y);
    mat(2,1) = gxy(x,y);
    mat(2,2) = gyy(x,y);
end

function [a,b,c] = getHessianABC(hessMat)
    a = hessMat(1,1);
    b = hessMat(1,2);
    c = hessMat(2,2);
end

function [eigvector1X, eigvector1Y, eigvector2X, eigvector2Y] = eigvectorsHessian(img)
    global Gx;
    global Gy;
    [xSize, ySize] = size(img);
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
    [xSize, ySize] = size(img);
    [iGxx, iGyy, iGxy] = getSecondDerivatives(img);
    hessGradientX = zeros(size(img),'like',img);
    hessGradientY = zeros(size(img),'like',img);
    for x = 1:xSize
        for y = 1:ySize
            localHessian = hessian(iGxx, iGyy, iGxy, x, y);
            [a,b,c] = getHessianABC(localHessian);
            hessGradientX(y, x) = a*Gx(y,x) + b*Gy(y,x);
            hessGradientY(y, x) = b*Gx(y,x) + c*Gy(y,x);
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
        isLin = determinant == 0;
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
    
    [xSize, ySize] = size(img);
    [iGxx, iGyy, iGxy] = getSecondDerivatives(img);
    vals = 0;
    for x = 1:xSize
        for y = 1:ySize
            determinant = Hx(y,x)*Gy(y,x) - Hy(y,x)*Gx(y,x);
       
            if(determinant == 0) && (~almostZero(Hx(y,x)) || ~almostZero(Hy(y,x))) && (~almostZero(Gx(y,x)) || ~almostZero(Gy(y,x)))
                hessVal = hessian(iGxx, iGyy, iGxy, x, y);
                eigVals = eig(hessVal);
                
                if ~(almostZero(eigVals(1)) || almostZero(eigVals(2)) || almostEqual(eigVals(1), eigVals(2)))
                    ptList(1,vals+1) = x;
                    ptList(2,vals+1) = y;
                    vals = vals+1;
                    %fprintf('determinant: %d, Hess: (%d, %d); Grad: (%d, %d), [%i,%i] \n', determinant, Hx(y,x), Hy(y,x), Gx(y,x), Gy(y,x), x, y);    
                end
                
            end
        end
    end
end

function ptList = heightPoints(img)
    global Gx;
    global Gy;
    [xSize, ySize] = size(img);
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
end

