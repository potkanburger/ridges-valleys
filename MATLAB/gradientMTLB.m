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

%old
%plotHeighPoints(ellipses_generated, ellipse_rotated, circle);

%last
plotRidgesValleys(ellipses_generated, ellipse_rotated, circle);
%plotFitLinRidgesValleys(ellipse_rotated);

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
    global Gx;
    global Gy;
    global Hx;
    global Hy;
    global EVs1;
    global EVs2;
    global closestEV1;
    global closestEV2;
    for i = 1:nargin   
        I = varargin{i}.data;
        
        HPG1 = getHeightPoints(I, 1);
        [ridges, pseudo_ridges, valleys, pseudo_valleys, excluded] = getRidgesValleys(HPG1, I, 10e-3);
        imgString = varargin{i}.name;
        
        figure(i)
        imshow(I);
        hold on;
        title(sprintf('Ridges and pseudo-ridges of %s', strrep(imgString, '_', ' ')));
        plotPoints(I, [ridges pseudo_ridges], 'ro');
        plotPoints(I, excluded, 'go');
        hold off;
        
        
        figure(nargin+i)
        imshow(I);
        hold on;
        title(sprintf('Valleys and pseudo-valleys of %s', strrep(imgString, '_', ' ')));
        plotPoints(I, [valleys pseudo_valleys], 'bo');
        plotPoints(I, excluded, 'go');
        saveas(gcf,strcat('valleys',imgString,'MATLAB.png'));
        hold off;
        
        
        clGx1 = onlyPx(excluded, Gx.*EVs1);
        clGy1 = onlyPx(excluded, Gy.*EVs1);
        clGx2 = onlyPx(excluded, Gx.*EVs2);
        clGy2 = onlyPx(excluded, Gy.*EVs2);
        clHx = onlyPx(excluded, Hx);
        clHy = onlyPx(excluded, Hy);
        
        
        figure(nargin*2+i)
        imshow(I);
        hold on;        
        title(sprintf('Quivers eigenvalues and plot of the chosen one for %s', strrep(imgString, '_', ' ')));
        quiver(clGx1, clGy1, 'b');
        quiver(clGx2, clGy2, 'r');
        quiver(clHx, clHy, 'g');
        plotPoints(I, closestEV1, 'bo');
        plotPoints(I, closestEV2, 'ro');
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

function plotFitLinRidgesValleys(varargin)
    for i = 1:nargin   
        I = varargin{i}.data;
        imgString = varargin{i}.name;
        [ridges, pseudo_ridges, valleys, pseudo_valleys, excluded] = getMathLinRidgesValleys(I, 9);
        
        figure(i)
        imshow(I);
        hold on;
        title(sprintf('Ridges and pseudo-ridges of %s', strrep(imgString, '_', ' ')));
        plotPoints(I, [ridges pseudo_ridges], 'ro');
        plotPoints(I, excluded, 'go');
        hold off;
        
        figure(nargin+i)
        imshow(I);
        hold on;
        title(sprintf('Valleys and pseudo-valleys of %s', strrep(imgString, '_', ' ')));
        plotPoints(I, [valleys pseudo_valleys], 'bo');
        plotPoints(I, excluded, 'go');
        hold off;
    end
end


function plotFitRidgesValleys(varargin)
    global Gx;
    global Gy;
    global Hx;
    global Hy;
    
    posV = nargin;
    
    for i = 1:nargin   
        I = varargin{i}.data;
        imgString = varargin{i}.name;
        [ridges, pseudo_ridges, valleys, pseudo_valleys, excluded] = getMathRidgesValleys(I, 12);
        
        figure(i)
        imshow(I);
        hold on;
        title(sprintf('Ridges and pseudo-ridges of %s', strrep(imgString, '_', ' ')));
        plotPoints(I, [ridges pseudo_ridges], 'ro');
        plotPoints(I, excluded, 'go');
        hold off;
        
        posV = posV+1;
        figure(posV)
        imshow(I);
        hold on;
        title(sprintf('Valleys and pseudo-valleys of %s', strrep(imgString, '_', ' ')));
        plotPoints(I, [valleys pseudo_valleys], 'bo');
        plotPoints(I, excluded, 'go');
        hold off;
    end
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
    
    global EVs1;
    global EVs2;
    
    determinants = zeros(size(img),'like',img);
    tolerances = zeros(size(img),'like',img);
    [ySize, xSize] = size(img);
    [iGxx, iGyy, iGxy] = getSecondDerivatives(img);
    
    %create EVs for tests
    EVs1 = zeros(size(img),'like',img);
    EVs2 = zeros(size(img),'like',img);
    for x = 1:xSize
        for y = 1:ySize
            hessVal = hessian(iGxx, iGyy, iGxy, x, y);
            eigVals = eig(hessVal);
            EVs1(y,x) = min(eigVals);
            EVs2(y,x) = max(eigVals);
        end
    end
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


function [rgList, pseudoRgList, valleyList, pseudoValleyList, excluded] = getRidgesValleys(heightPoints,img, d_tolerance)
    global Gx;
    global Gy;
    global Hx;
    global Hy;
    
    global closestEV1;
    global closestEV2;
    
    closestEV1 = [];
    closestEV2 = [];
    valsCSE1 = 0;
    valsCSE2 = 0;
    
    nbHp = size(heightPoints);
    nbHp = nbHp(2);
    [iGxx, iGyy, iGxy] = getSecondDerivatives(img);
    valsRg = 0;
    valsPsRg = 0;
    valsValley = 0;
    valsPsValley = 0;
    valsExc = 0;
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
        
        %lengthH = sqrt(HxVal*HxVal + HyVal*HyVal);
        %lengthG = sqrt(GxVal*GxVal + GyVal*GyVal);
        
        diffEV1 = abs(HxVal - ev1*GxVal) + abs(HyVal - ev1*GyVal);
        diffEV2 = abs(HxVal - ev2*GxVal) + abs(HyVal - ev2*GyVal);
        %diffEV1 = abs(lengthH - ev1*lengthG);
        %diffEV2 = abs(lengthH - ev2*lengthG);
        %idea?
        isEqEV1 = diffEV1 < diffEV2 && diffEV1 < d_tolerance;
        isEqEV2 = diffEV2 < diffEV1 && diffEV2 < d_tolerance;
        
        %isEqEV1 = diffEV1 < diffEV2;
        %isEqEV2 = diffEV2 < diffEV1;
        
        if(diffEV1 < diffEV2)
            closestEV1(1,valsCSE1+1) = xVal;
            closestEV1(2,valsCSE1+1) = yVal;
            valsCSE1 = valsCSE1+1;
        elseif(diffEV2 < diffEV1)
            closestEV2(1,valsCSE2+1) = xVal;
            closestEV2(2,valsCSE2+1) = yVal;
            valsCSE2 = valsCSE2+1;
        end
        
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
        else
            excluded(1,valsExc+1) = xVal;
            excluded(2,valsExc+1) = yVal;
            valsExc = valsExc+1;
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
    
    if(valsExc == 0)
        excluded = [];
    end
    
    if(valsCSE1 == 0)
        closestEV1 = [];
    end
    
    if(valsCSE2 == 0)
        closestEV2 = [];
    end
end


function [rgList, pseudoRgList, valleyList, pseudoValleyList, undetermined] = getMathRidgesValleys(img, area_sidesize)
    valsRg = 0;
    valsPsRg = 0;
    valsValley = 0;
    valsPsValley = 0;
    valsExc = 0;


    if mod(area_sidesize,2) == 0
        area_sidesize = area_sidesize + 1;
    end
    
    area_half = (area_sidesize-1)/2; 
    
    
    calcImg = img;
    f_col = calcImg(:, 1);
    l_col = calcImg(:,end);
    for i = 1:area_half
        calcImg = [f_col calcImg l_col];
    end
    
    f_line = calcImg(1, :);
    l_line = calcImg(end,:);
    
    for i = 1:area_half
        calcImg = [f_line; calcImg;l_line];
    end
    
    
    
    local_area_data = zeros([area_sidesize, area_sidesize], 'like', img); 
       
    [ySize, xSize] = size(calcImg);
    
    startPt = 1 + area_half;
    
    xValues = [];
    for t = 1:area_sidesize
       xValues = [xValues repmat(t, 1, area_sidesize)];
    end
    xValues = [xValues]';
    
    yValues = repmat([1:area_sidesize],1,area_sidesize);
    yValues = [yValues]';
    
    nbValsArray = size(yValues);
    nbVals = nbValsArray(1);
    
    localCenter = 1+area_half;
    opts = optimoptions(@fsolve,'Display','off', 'Algorithm', 'levenberg-marquardt', 'TolFun', 1e-12, 'TolX', 1e-12);
    
    for x = startPt:(xSize-area_half)
        for y = startPt:(ySize-area_half)
            
            for localX = 1:area_sidesize
                imgX = x - area_half + localX - 1;
                for localY = 1:area_sidesize
                    imgY = y - area_half + localY - 1;
                    local_area_data(localY,localX) = calcImg(imgY, imgX);
                end
            end
            
            
            zValues = zeros(nbVals,1);
            for tmp = 1:nbVals
                zValues(tmp) = local_area_data(yValues(tmp), xValues(tmp));
            end
            
            
            
            sf = fit([xValues, yValues],zValues,'poly23');
            coeffs = coeffvalues(sf);
            p00 = coeffs(1);
            p10 = coeffs(2);
            p01 = coeffs(3);
            p20 = coeffs(4);
            p11 = coeffs(5);
            p02 = coeffs(6);
            p21 = coeffs(7);
            p12 = coeffs(8);
            p03 = coeffs(9);
            
            syms f(a,b);
            f(a,b) = p00 + p10*a + p01*b + p20*a^2 + p11*a*b + p02*b^2 + p21*a^2*b + p12*a*b^2 + p03*b^3;
            
            gx = diff(f, a);
            gy = diff(f, b);
            
            gxx = diff(gx, a);
            gxy = diff(gx, b);
            
            gyx = diff(gy, a);
            gyy = diff(gy, b);
            
            Hgx = gxx*gx + gxy*gy;
            Hgy = gyx*gx + gyy*gy;
            
            determinant = Hgx*gy - Hgy*gx;
            eqn = @(p) double(determinant(p(1),p(2)));
            
            signNb = sign(double(determinant(localCenter-1, localCenter-1)));
            changeSign = false;
            for sX = 0:2
                for sY = 0:2
                    if(sign(double(determinant(localCenter-1+sY, localCenter-1+sX))) ~= signNb)
                        changeSign = true;                        
                    end
                end
            end
            
            if(changeSign)
                S = fsolve(eqn, [localCenter-1,localCenter-1], opts);


                Hessian = [gxx(S(1), S(2)), gxy(S(1), S(2)); gyx(S(1), S(2)), gyy(S(1), S(2))];
                Hg = [Hgx(S(1), S(2)), Hgy(S(1), S(2))];
                Grad = [gx(S(1), S(2)), gy(S(1), S(2))];

                eigVals = eig(Hessian);

                ev1 = min(eigVals);
                ev2 = max(eigVals);

                diffEV1 = abs(Hg(1)-Grad(1)*ev1) + abs(Hg(2)-Grad(2)*ev1);
                diffEV2 = abs(Hg(1)-Grad(1)*ev2) + abs(Hg(2)-Grad(2)*ev2);

                xVal = round(S(1));
                yVal = round(S(2));
                
                xInCalcImg = x-area_half + xVal-1;
                yInCalcImg = y-area_half + yVal-1;
                
                xInInputImg = xInCalcImg -area_half;
                yInInputImg = yInCalcImg -area_half;
                
                if(xInInputImg >= 1 && xInInputImg <= xSize - 2*area_half && yInInputImg >= 1 && yInInputImg <= ySize - 2*area_half)
                
                    if(diffEV1 < diffEV2 && diffEV1 <10e-6)
                        if ev2 > 0
                            %valley
                            valleyList(1,valsValley+1) = xInInputImg;
                            valleyList(2,valsValley+1) = yInInputImg;
                            valsValley = valsValley+1;
                        elseif ev2 < 0
                            %pseudo-ridge
                            pseudoRgList(1,valsPsRg+1) = xInInputImg;
                            pseudoRgList(2,valsPsRg+1) = yInInputImg;
                            valsPsRg = valsPsRg+1;
                        end
                    elseif(diffEV2 < diffEV1 && diffEV2 <10e-6)
                        if ev1 < 0
                            %ridge
                            rgList(1,valsRg+1) = xInInputImg;
                            rgList(2,valsRg+1) = yInInputImg;
                            valsRg = valsRg+1;
                        elseif ev1 > 0
                            %pseudo-valley
                            pseudoValleyList(1,valsPsValley+1) = xInInputImg;
                            pseudoValleyList(2,valsPsValley+1) = yInInputImg;
                            valsPsValley = valsPsValley+1;
                        end
                    else                
                    undetermined(1,valsExc+1) = xInInputImg;
                    undetermined(2,valsExc+1) = yInInputImg;
                    valsExc = valsExc+1;
                    end
                end
            end
%             [fx, fy, fxx, fxy, fyy] = differentiate( sf, localCenter, localCenter);
%             g = [fx, fy];
%             H = [fxx, fxy; fxy, fyy];
%             eigVals = eig(H);
%             
%             ev1 = min(eigVals);
%             ev2 = max(eigVals);
%             
%             t_g = g.';
%             Hg = H * t_g;
%             
%             determinant = Hg(1)*g(2) - Hg(2)*g(1)
%             diff_valEV1x = abs(Hg(1)-ev1*g(1));
%             diff_valEV1y = abs(Hg(2)-ev1*g(2));
%             
%             diff_valEV2x = abs(Hg(1)-ev2*g(1));
%             diff_valEV2y = abs(Hg(2)-ev2*g(2)); 
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
    
    if(valsExc == 0)
        undetermined = [];
    end
    
end


function [rgList, pseudoRgList, valleyList, pseudoValleyList, undetermined] = getMathLinRidgesValleys(img, area_sidesize)
    valsRg = 0;
    valsPsRg = 0;
    valsValley = 0;
    valsPsValley = 0;
    valsExc = 0;
    
    if mod(area_sidesize,2) == 0
        area_sidesize = area_sidesize + 1;
    end
    
    area_half = (area_sidesize-1)/2; 
    
    
    calcImg = img;
    f_col = calcImg(:, 1);
    l_col = calcImg(:,end);
    for i = 1:area_half
        calcImg = [f_col calcImg l_col];
    end
    
    f_line = calcImg(1, :);
    l_line = calcImg(end,:);
    
    for i = 1:area_half
        calcImg = [f_line; calcImg;l_line];
    end
    
    [ySize, xSize] = size(calcImg);
    
    local_area_data_x = zeros([1, area_sidesize], 'like', img)'; 
    local_area_data_y = zeros([1, area_sidesize], 'like', img)';
    
    xyValues = [1:area_sidesize]';
    
    localCenter = 1+area_half;
    nbVals = area_sidesize;
    startPt = 1 + area_half;
    
    
    for x = startPt:(xSize-area_half)
        for y = startPt:(ySize-area_half)
        
            for localX = 1:area_sidesize
                imgX = x - area_half + localX - 1;
                local_area_data_x(localX) = calcImg(y, imgX);
            end
            
            for localY = 1:area_sidesize
                imgY = y - area_half + localY - 1;
                local_area_data_y(localY) = calcImg(imgY, x);
            end
            
            xFuncFit = fit(xyValues,local_area_data_x,'poly5');
            yFuncFit = fit(xyValues,local_area_data_y,'poly5');
            
            syms a b;
          
            coeffsX = coeffvalues(xFuncFit);
            p1 = coeffsX(1);
            p2 = coeffsX(2);
            p3 = coeffsX(3);
            p4 = coeffsX(4);
            p5 = coeffsX(5);
            p6 = coeffsX(6);
            xFunc = p1*a^5 + p2*a^4 + p3*a^3 + p4*a^2 + p5*a + p6;
            
            coeffsY = coeffvalues(yFuncFit);
            p1 = coeffsY(1);
            p2 = coeffsY(2);
            p3 = coeffsY(3);
            p4 = coeffsY(4);
            p5 = coeffsY(5);
            p6 = coeffsY(6);
            
            yFunc = p1*b^5 + p2*b^4 + p3*b^3 + p4*b^2 + p5*b + p6;
            
            gx = diff(xFunc);
            gxx = diff(xFunc,2);
            
            gy = diff(yFunc);
            gyy = diff(yFunc,2);
            
            Hgx = gxx*gx;
            Hgy = gyy*gy;
            
            determinantX = @(p) double(subs(Hgx, a, p));
            determinantY = @(p) double(subs(Hgy, b, p));
            deter = Hgx*gy - Hgy*gx;
            deterX = @(p) subs(deter, a, p);
            determinant = @(m,n) double(subs(deterX(n), b, m));

            signNb = sign(determinant(localCenter-1, localCenter-1));
            changeSign = false;
            for sX = 0:2
                for sY = 0:2
                    if(sign(determinant(localCenter-1+sY, localCenter-1+sX)) ~= signNb)
                        changeSign = true;                        
                    end
                end
            end
            
            
            if(changeSign)
 
                xVal = localCenter;
                yVal = localCenter;
                
                xInCalcImg = x-area_half + xVal-1;
                yInCalcImg = y-area_half + yVal-1;
                
                xInInputImg = xInCalcImg -area_half;
                yInInputImg = yInCalcImg -area_half;
                
                undetermined(1,valsExc+1) = xInInputImg;
                undetermined(2,valsExc+1) = yInInputImg;
                valsExc = valsExc+1;
            
            else
                res = sprintf("x,y: %d,%d", x, y);
                
                for sX = 0:2
                    for sY = 0:2
                        res = [res, ',', sprintf("%d", determinant(localCenter-1+sY,localCenter-1+sX))];
                    end
                end
                
                res= join(res)
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
    
    if(valsExc == 0)
        undetermined = [];
    end
    
end

