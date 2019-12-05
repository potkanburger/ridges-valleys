
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
inverted_ellipse_rotated = invertImg(ellipse_rotated);

ITEMP = double(rr);
circle = createImgStruct('circle', ITEMP, true);

ITEMP = simpleComponentTest();
cpTest = createImgStruct('componentTest', ITEMP, true);

global ParTree;
global RnkTree;
global ParNode;
global RnkNode;

global nodes;
global lowestNode;

[r, nd, cpMap] = makeComponentTree(cpTest);

lM = getLocalMaxima(r, nd);
nbLocalMaxima = size(lM(:), 1);
figureCount = 1;

figure(figureCount)
hold on;
X = mat2gray(cpTest.data, [0 255]);
imshow(X, 'InitialMagnification','fit');
hold off;

for max = 1:nbLocalMaxima
   path{max} = getComponentPath(lM(max), r, nd);
   
   curPath = path{max};
   pixels = [];
   pixels(:) = uint16(0);
   for p = 1:size(curPath(:),1)
        figureCount = figureCount + 1;
        cpId = curPath(p);
        pixels = union(pixels, getPixelList(cpId, cpMap));
        lvlImg = [];
        lvlImg(cpTest.ySize, cpTest.xSize) = 0;
        for px = 1:size(pixels(:), 1)
            [y, x] = getYXfromID(cpTest, pixels(px));
            lvlImg(y, x) = 1;
        end
       
        figure(figureCount)
        hold on;
        X = mat2gray(lvlImg);
        imshow(X, 'InitialMagnification','fit');
        hold off;
       
   end
    
   
end




% firstLevelParents = [];
% firstLevelParents(:) = uint16(0);
% for max = 1:size(lM, 1)
%    parents = getParents(nd, lM(max));
%    firstLevelParents = union(firstLevelParents, parents);
% end
% 
% pixels = [];
% pixels(:) = uint16(0);
% for cp = 1:size(firstLevelParents(:), 1)
%     cpId = firstLevelParents(cp);
%     pixels = union(pixels, getPixelList(cpId, cpMap));
% end
% pixels = union(pixels, lM);
% 
% figure(1)
% hold on;
% X = mat2gray(cpTest.data, [0 255]);
% imshow(X, 'InitialMagnification','fit');
% hold off;
% 
% localMaxImg(cpTest.ySize, cpTest.xSize) = 0;
% for max = 1:size(lM, 1)
%     [y, x] = getYXfromID(cpTest, lM(max));
%     localMaxImg(y, x) = 1;
% end


%functions

function [root, nodesRes, componentMap] = makeComponentTree(IMG)
    global ParTree;
    global RnkTree;
    global ParNode;
    global RnkNode;

    global nodes;
    global lowestNode;
    
    MakeSets(IMG);
    initNodes(IMG);
    pxOrder = getPixelOrder(IMG);
    pxSize = size(pxOrder, 1);

    componentMap = uint16(1:pxSize);

    for pxNb = 1:pxSize
        curId = pxOrder(pxNb);
        curTree = FindInTree(curId);
        curNode = FindInNode(lowestNode(curTree));
        listValidNgb = getAlreadyProcessedNeighbors(IMG, pxOrder, pxNb);
        nbNgb = size(listValidNgb,2);
        if(nbNgb > 0)
            for ngbLocalId = 1:nbNgb
                nId = listValidNgb(ngbLocalId);
                adjTree = FindInTree(nId);
                adjNode = FindInNode(lowestNode(adjTree));

                if(curNode ~= adjNode)
                    if(nodes(curNode).Level == nodes(adjNode).Level)
                        curNode = MergeNodes(adjNode, curNode);
                    else
                        nodes(curNode).Children = union(nodes(curNode).Children, adjNode);
                        nodes(curNode).Area = nodes(curNode).Area + nodes(adjNode).Area;
                        nodes(curNode).Highest = max(nodes(curNode).Highest, nodes(adjNode).Highest);
                    end
                    curTree = LinkTree(adjTree, curTree);
                    lowestNode(curTree) = curNode;
                end
            end
        end
    end
    root = lowestNode(FindInTree(FindInNode(1)));

    for pxNb = 1:pxSize
        componentMap(pxNb) = FindInNode(pxNb);
    end
    
    nodesRes = nodes;
end

function MakeSets(img)
    global ParTree;
    global RnkTree;
    global ParNode;
    global RnkNode;
    
    nbPx = img.xSize*img.ySize;
    RnkNode = uint16(0);
    RnkTree = uint16(0);    
    RnkNode(nbPx) = 0;
    RnkTree(nbPx) = 0;
    ParTree = uint16(1:nbPx);
    ParNode = uint16(1:nbPx);
end

function initNodes(img)
    global nodes;
    global lowestNode;
    
    nbPx = img.xSize*img.ySize;
    nodes = cpTreeNode();
    nodes(nbPx) = cpTreeNode();
    lowestNode = uint16(1:nbPx);
    for px = 1:nbPx
        nodes(px) = MakeNode(getVal(img, px));
    end
    
end

function node = MakeNode(level)
    node = cpTreeNode();
    node.Level = level;
    node.Highest = level;
end

function resNode = MergeNodes(node1, node2)
    global nodes;
    resNode = LinkNode(node1, node2);
    if(resNode == node2)
        nodes(node2).Children = union(nodes(node2).Children, nodes(node1).Children);
    else
        nodes(node1).Children = union(nodes(node2).Children, nodes(node1).Children);
    end
    
    nodes(resNode).Area = nodes(node2).Area + nodes(node1).Area;
    nodes(resNode).Highest = max(nodes(node1).Highest, nodes(node2).Highest);    
end


function index = LinkTree(elemX, elemY)
    global RnkTree;
    global ParTree;
    x = elemX;
    y = elemY;
    rankX = RnkTree(x);
    rankY = RnkTree(y);
    
    if(rankX > rankY)
        x = elemY;
        y = elemX;
    elseif(rankX == rankY)
        RnkTree(y) = rankY + 1;
    end
    
    ParTree(x) = y;    
    index = y;    
end

function index = LinkNode(elemX, elemY)
    global RnkNode;
    global ParNode;
    x = elemX;
    y = elemY;
    rankX = RnkNode(x);
    rankY = RnkNode(y);
    
    if(rankX > rankY)
        x = elemY;
        y = elemX;
    elseif(rankX == rankY)
        RnkNode(y) = rankY + 1;
    end
    
    ParNode(x) = y;    
    index = y; 
end

function index = FindInTree(id)
    global ParTree;
    parId = ParTree(id);
    if(parId ~= id)
        ParTree(id) = FindInTree(parId);
    end
    index = ParTree(id);
end

function index = FindInNode(id)
    global ParNode;
    parId = ParNode(id);
    if(parId ~= id)
        ParNode(id) = FindInNode(parId);
    end
    index = ParNode(id);
end

function index = getIndex(img, x, y)
    index = img.xSize*(y-1)+ x;
end

function [y, x] = getYXfromID(img, id)
    y = ceil(double(id)/img.xSize);
    x = id - img.xSize*(y-1);
end

function pixels = getPixelList(componentId, cpMap)
    nbPxMax = size(cpMap, 2);
    pixels(nbPxMax) = uint16(0);
    nbPx = 0;
    for px = 1:nbPxMax
        if(cpMap(px) == componentId)
            nbPx = nbPx + 1;
            pixels(nbPx) = px;
        end
    end
    
    if(nbPx == 0)
        pixels = [];
        pixels(:) = uint16.empty;
    else
        pixels = pixels(1:nbPx);
    end    
end

function value = getVal(img, id)
    [yVal, xVal] = getYXfromID(img,id);
    value = img.data(yVal,xVal);
end

function orderedPixels = getPixelOrder(img)
    nbPx = img.xSize*img.ySize;
    tmpList(nbPx, :) = [0 0]; 
    for px = 1:nbPx
        tmpList(px,:) = [px getVal(img, px)];
    end
    result = sortrows(tmpList, 2, 'descend');
    orderedPixels = result(:,1);        
end

function validNbList = getAlreadyProcessedNeighbors(img, orderedPixels, curPos)
    %4-nbhood, valid -> no test for value because if already processed, the
    %value is >= cur
    nbs = 0;
    id = orderedPixels(curPos);
    processed = orderedPixels(1:curPos);
    
    [yVal, xVal] = getYXfromID(img,id);
    if(xVal > 1)
        testIndex = id - 1;
        if(ismember(testIndex, processed))
            nbs = nbs+1;
            validNbList(nbs) = testIndex;
        end
    end
    
    if(xVal < img.xSize)
        testIndex = id + 1;
        if(ismember(testIndex, processed))
            nbs = nbs+1;
            validNbList(nbs) = testIndex;
        end
    end
    
    if(yVal > 1)
        testIndex = id - img.xSize;
        if(ismember(testIndex, processed))
            nbs = nbs+1;
            validNbList(nbs) = testIndex;
        end
    end
    
    
    if(yVal < img.ySize)
        testIndex = id + img.xSize;
        if(ismember(testIndex, processed))
            nbs = nbs+1;
            validNbList(nbs) = testIndex;
        end
    end
    
    if(nbs == 0)
        validNbList = {};
    end
    
end

function children = getChildren(node, nodes)
    children = nodes(node).Children;
end

function localMax = getLocalMaxStep(nd, nodes)
    c = nodes(nd).Children;
    childrenCount = size(c,1);
    
    localMax = [];
    localMax(:) = uint16.empty;
    
    if(childrenCount == 0)
        localMax = union(localMax, nd);
    else
        for cd = 1:childrenCount
            tmpLocalMax = getLocalMaxStep(c(cd), nodes);
            localMax = union(localMax, tmpLocalMax);           
        end
    end
end

function localMaxList = getLocalMaxima(root, nodes)
    localMaxList = getLocalMaxStep(root, nodes);
end

function path = getLeaf(target, element, nodes, curPath)
    if(element == target)
        path = [curPath element];
    else
        elements = nodes(element).Children;
        elemSize = size(elements(:), 1);
        if(elemSize == 0)
            path = [];
        else
            localPath = [curPath element];
            path = [];
            for e = 1:elemSize
                testPath =  getLeaf(target, elements(e), nodes, localPath);
                if(size(testPath(:), 1) > 0 && ismember(target, testPath))
                    path = testPath;
                    break;
                end
            end
        end
    
    end
end


%cp path from pt to root
function path = getComponentPath(pt, root, nodes)
    curPath = [];
    curPath(:) = uint16.empty;    
    path = flip(getLeaf(pt, root, nodes, curPath));
end

function parents = getParents(nodes, node)
    parCount = 0;
    sizeNodes = size(nodes,2);
    parents(sizeNodes) = uint16(0);
    
    for pt = 1:sizeNodes
        if(ismember(node, nodes(pt).Children))
            parCount = parCount+1;
            parents(parCount) = pt;
        end
    end
    
    if(parCount == 0)
        parents = [];
        parents(:) = uint16.empty;
    else
        parents = parents(1:parCount);
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


function res = simpleComponentTest()
    sizeX = 10;
    sizeY = 8;
    res(sizeY, sizeX) = 0;
    res(1:sizeY-5, 1:sizeX-5) = 10;
    res(sizeY-4:sizeY, :) = 50;
    res(:, sizeX-4:sizeX) = 50;
    res(sizeY-3, sizeX-3) = 100;
    
    res(sizeY-4, sizeX-4) = 80;
    res(sizeY-4, sizeX-3) = 80;
    res(sizeY-4, sizeX-2) = 80;
    
    res(sizeY-3, sizeX-4) = 80;
    res(sizeY-3, sizeX-2) = 80;
    
    res(sizeY-2, sizeX-4) = 80;
    res(sizeY-2, sizeX-3) = 80;
    res(sizeY-2, sizeX-2) = 80;
    
    res(sizeY,sizeX) = 130;
end


function inverted = invertImg(imgStruct)
    inverted = imgStruct;
    inverted.name = "inverted_"+imgStruct.name;
    inverted.data = -imgStruct.data;
end

function img = createImgStruct(name, data, isGray)
    img.name = name;
    if(~isGray)
        img.data = rgb2gray(data);
    else
        img.data = data;
    end
    [ySize, xSize] = size(img.data);
    img.xSize = xSize;
    img.ySize = ySize;
    
end