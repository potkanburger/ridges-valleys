imgString = 'ellipses_generated';
imgExtension = '.jpg';
fullImgPath = strcat(imgString,imgExtension);

ITEMP = imread(fullImgPath);
I = ITEMP;
% I = rgb2gray(ITEMP);

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
close all

function res = almostZero(val)
    res = abs(val) <= 0.00000000001;        
end

function res = almostEqual(val1, val2)
    res = almostZero(val1-val2);     
end

function plotPoints(img, pts)
    xValues = pts(1,:);
    yValues = pts(2,:);
    plot(xValues, yValues, 'r+', 'LineWidth', 2, 'MarkerSize', 2);
end

function [gxx, gyy, gxy] = getSecondDerivatives(img)
  [gx, gy] = gradient(double(img));
  [gxx, gxy] = gradient(double(gx));
  [gyx, gyy] = gradient(double(gy));
end

function mat = hessian(gxx, gyy, gxy, x, y)
    mat(1,1) = gxx(x,y);
    mat(1,2) = gxy(x,y);
    mat(2,1) = gxy(x,y);
    mat(2,2) = gyy(x,y);
end

function [a,b,c] = getHessianABC(hessMat)
    a = hessMat(1,1)
    b = hessMat(1,2)
    c = hessMat(2,2)
end

function isLin = areLambdaLinear(hessMat, gxVal, gyVal)
    [a,b,c] = getHessianABC(hessMat);
    if ~almostZero(gyVal) && ~almostZero(gxVal)
        lambda1 = (a*gyVal + b*gxVal)/gyVal
        lambda2 = (b*gyVal + c*gxVal)/gxVal
        
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
    [iGxx, iGyy, iGxy] = getSecondDerivatives(img)
    for x = 1:xSize
        for y = 1:ySize
            gxValue = iGx(x,y);
            gyValue = iGy(x,y);       
            if almostZero(gxValue) || almostZero(gyValue)
                hessVal = hessian(iGxx, iGyy, iGxy, x, y);
                eigVals = eig(hessVal);
                if ~(almostZero(eigVals(1)) || almostZero(eigVals(2)) || almostEqual(eigVals(1), eigVals(2)))
                    if areLambdaLinear(hessVal, gxValue, gyValue)
                        ptList(1,vals+1) = x;
                        ptList(2,vals+1) = y;
                        vals = vals+1;
                    end
                end
            end
        end
    end
end

