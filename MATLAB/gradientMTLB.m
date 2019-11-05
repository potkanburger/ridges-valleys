% parameters
clear all;


imgString = 'ellipses_generated';
imgExtension = '.jpg';
gray = true;

%program
fullImgPath = strcat(imgString,imgExtension);
ITEMP = imread(fullImgPath);
ITEMP = im2double(ITEMP);

ellipses_generated = createImgStruct(imgString, ITEMP, true);


ITEMP = generateEllipse(50,60);

ellipse_rotated = createImgStruct('rotated_ellipse', ITEMP, true);

ITEMP = double(rr);
circle = createImgStruct('circle', ITEMP, true);

%plotHeighPoints(ellipses_generated, ellipse_rotated, circle);
plotRidgesValleys(ellipses_generated, ellipse_rotated, circle);


%functions

function img = createImgStruct(name, data, isGray)
    img.name = name;
    if(~isGray)
        img.data = rgb2gray(data);
    else
        img.data = data;
    end
end

function plotHeighPoints(varargin)
    for i = 1:nargin
        figure(i)
        I = varargin{i}.data;
        imshow(I);
        hold on;
        HPG1 = getHeightPoints(I, 1);
        imgString = varargin{i}.name;
        title(sprintf('Height Points of %s', strrep(imgString, '_', ' ')));
        plotPoints(I, HPG1, 'ro');
        saveas(gcf,strcat('hpts_',imgString,'MATLAB.png'));
        hold off;
    end
end

function plotRidgesValleys(varargin)
    for i = 1:nargin   
        I = varargin{i}.data;
        HPG1 = getHeightPoints(I, 1);
        [ridges, pseudo_ridges, valleys, pseudo_valleys] = getRidgesValleys(HPG1, I);
        imgString = varargin{i}.name;
        
        figure((i-1)*2 + 1)
        imshow(I);
        hold on;
        title(sprintf('Ridges and pseudo-ridges of %s', strrep(imgString, '_', ' ')));
        plotPoints(I, [ridges pseudo_ridges], 'ro');
        saveas(gcf,strcat('ridges_',imgString,'MATLAB.png'));
        hold off;
        
        figure(i*2)
        imshow(I);
        hold on;
        title(sprintf('Valleys and pseudo-valleys of %s', strrep(imgString, '_', ' ')));
        plotPoints(I, [valleys pseudo_valleys], 'ro');
        saveas(gcf,strcat('valleys',imgString,'MATLAB.png'));
        hold off;
        
        
        
    end
end

function res = generateEllipse(sizeX, sizeY)

alpha = pi*30/180;
a = 20;
b = 30;
sx = 0;
sy = 5;

ellipse = sqrt(((xx(sizeX,sizeY) - sx)*cos(alpha) + (yy(sizeX,sizeY) -sy)*sin(alpha))^2/a^2 + ((xx(sizeX,sizeY) - sx)*sin(alpha) - (yy(sizeX,sizeY) -sy)*cos(alpha))^2/b^2);

res = double(ellipse);
end


function ptList = getHeightPoints(img, sigma)
    global Gx;
    global Gy;
    global Hx;
    global Hy;

    [Gx,Gy] = gaussfiltgradientxy(img,sigma);
    [Hx, Hy] = hessianXgradient(img);
    
    ptList = calculateHeightPoints(img);
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
            tt = [Gx(y,x) Gy(y,x)];
            tt = tt.';
            test = localHessian * tt;
            
            hessGradientX(y, x) = test(1);
            hessGradientY(y, x) = test(2);
        end
    end
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

function ptList = calculateHeightPoints(img)
    global Gx;
    global Gy;
    
    global Hx;
    global Hy;
    
    determinants = zeros(size(img),'like',img);
    tolerances = zeros(size(img),'like',img);
    [ySize, xSize] = size(img);
    [iGxx, iGyy, iGxy] = getSecondDerivatives(img);
    vals = 0;
    for x = 1:xSize
        previousSign = 0;
        previousDet = 0;        
        for y = 1:ySize
            [determinant, tol] = getDeterminant(y,x);
            determinants(y, x) = determinant;
            tolerances(y,x) = tol;
            diffEV1 = 0;
            diffEV2 = 0;
            
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
            determinant = determinants(y,x);
            tol = tolerances(y,x);
            
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


function [rgList, pseudoRgList, valleyList, pseudoValleyList] = getRidgesValleys(heightPoints, img)
    global Gx;
    global Gy;
    global Hx;
    global Hy;
    nbHp = size(heightPoints);
    nbHp = nbHp(2);
    [iGxx, iGyy, iGxy] = getSecondDerivatives(img);
    valsRg = 0;
    valsPsRg = 0;
    valsValley = 0;
    valsPsValley = 0;
    for pt = 1:nbHp
        xVal = heightPoints(1, pt);
        yVal = heightPoints(2, pt);
        
        hessVal = hessian(iGxx, iGyy, iGxy, xVal, yVal);
        eigVals = eig(hessVal);
        
        HxVal = Hx(yVal,xVal);
        HyVal = Hy(yVal,xVal);
        GxVal = Gx(yVal,xVal);
        GyVal = Gy(yVal,xVal);
        
        ev1 = min(eigVals);
        ev2 = max(eigVals);
        
        lengthH = sqrt(HxVal*HxVal + HyVal*HyVal);
        lengthG = sqrt(GxVal*GxVal + GyVal*GyVal);

        diffEV1 = abs(lengthH - ev1*lengthG);
        diffEV2 = abs(lengthH - ev2*lengthG);
        %idea?
        %isEqEV1 = diffEV1 < diffEV2 && diffEV1 < lengthH/1000;
        %isEqEV2 = diffEV2 < diffEV1 && diffEV2 < lengthH/1000;
        
        isEqEV1 = diffEV1 < diffEV2;
        isEqEV2 = diffEV2 < diffEV1;
        
        if isEqEV1
            if ev2 > 0
                %valley
                valleyList(1,valsValley+1) = xVal;
                valleyList(2,valsValley+1) = yVal;
                valsValley = valsValley+1;
            elseif ev2 < 0
               %pseudo-ridge
                pseudoRgList(1,valsPsRg+1) = xVal;
                pseudoRgList(2,valsPsRg+1) = yVal;
                valsPsRg = valsPsRg+1;
            end
        elseif isEqEV2
            if ev1 < 0
                %ridge
                rgList(1,valsRg+1) = xVal;
                rgList(2,valsRg+1) = yVal;
                valsRg = valsRg+1;
            elseif ev1 > 0
                %pseudo-valley
                pseudoValleyList(1,valsPsValley+1) = xVal;
                pseudoValleyList(2,valsPsValley+1) = yVal;
                valsPsValley = valsPsValley+1;
            end
        end
    end
    if(valsRg == 0)
        rgList = [];
    end
    
    if(valsValley == 0)
        valleyList = [];
    end
    
    if(valsPsRg == 0)
        pseudoRgList = [];
    end
    
    if(valsPsValley == 0)
        pseudoValleyList = [];
    end
end
