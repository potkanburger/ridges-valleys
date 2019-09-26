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

[Gx,Gy] = imgradientxy(I, 'central');
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
plotPoints(I, mL);
title('HeightPoints');
saveas(gcf,strcat('hpts_',imgString,'MATLAB.png'))
hold off

[Hx, Hy] = hessianXgradient(I);

figure(4)
imshow(I);
hold on;
plotPointsBig(I, mL);
title('Gradient Vectors');
quiver(Gx, Gy);
quiver(Hx, Hy);

saveas(gcf,strcat('vectors_',imgString,'_MATLAB.png'))
hold off
% close all

function res = almostZero(val)
    res = abs(val) <= 0.000000000001;   
end

function res = almostEqual(val1, val2)
    res = almostZero(val1-val2);     
end

function plotPoints(img, pts)
    xValues = pts(1,:);
    yValues = pts(2,:);
    plot(xValues, yValues, 'r+', 'LineWidth', 2, 'MarkerSize', 2);
end

function plotPointsBig(img, pts)
    xValues = pts(1,:);
    yValues = pts(2,:);
    plot(xValues, yValues, 'ro', 'LineWidth', 4, 'MarkerSize', 4);
end


function [gxx, gyy, gxy] = getSecondDerivatives(img)
  [gx, gy] = imgradientxy(double(img), 'central');
  [gxx, gxy] = imgradientxy(double(gx), 'central');
  [gyx, gyy] = imgradientxy(double(gy), 'central');
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

function [hessGradientX,hessGradientY] = hessianXgradient(img)
    [xSize, ySize] = size(img);
    [iGx,iGy] = imgradientxy(img, 'central');
    [iGxx, iGyy, iGxy] = getSecondDerivatives(img);
    hessGradientX = zeros(size(img),'like',img);
    hessGradientY = zeros(size(img),'like',img);
    for x = 1:xSize
        for y = 1:ySize
            localHessian = hessian(iGxx, iGyy, iGxy, x, y);
            [a,b,c] = getHessianABC(localHessian);
            hessGradientX(x, y) = a*iGx(x,y) + b*iGy(x,y);
            hessGradientY(x, y) = b*iGx(x,y) + c*iGy(x,y);
        end
    end
end

function isLin = pixelIsLinear(hessVal, gxValue, gyValue)
    [a,b,c] = getHessianABC(hessVal);
    hessGradientX = a*gxValue + b*gyValue;
    hessGradientY = b*gxValue + c*gyValue;
    
    determinant = hessGradientX*gyValue - hessGradientY*gxValue;    
    isLin = determinant == 0 && ~(almostZero(hessGradientX*gyValue) && almostZero(hessGradientY*gxValue));
    
    if isLin
        hessGradientX
        hessGradientY
        gxValue
        gyValue
    end
    
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

function ptList = heightPoints(img)
    [xSize, ySize] = size(img);
    [iGx,iGy] = imgradientxy(img, 'central');
    vals = 0;
    [iGxx, iGyy, iGxy] = getSecondDerivatives(img);
    for x = 1:xSize
        for y = 1:ySize
            gxValue = iGx(x,y);
            gyValue = iGy(x,y);       
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

