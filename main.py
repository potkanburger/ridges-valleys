import cv2 as cv
import numpy as np
import math
#import matplotlib.pyplot as plot

def getImgGrayscale(url):
    return cv.imread(url, 0).astype(np.float)

def createImg(width, height):
    return np.zeros((height, width, 1), dtype="float")

def sqr(val):
    return pow(val, 2)

def xor(bool1, bool2):
    return bool1 ^ bool2


def addNeumannBorder(img, size):
    return cv.copyMakeBorder(img, size, size, size, size, cv.BORDER_REPLICATE)



def hessian(img):
    width = np.size(img, 0)
    height = np.size(img, 1)
    hessianM = {}
    for x in range(1, width - 2):
        for y in range(1, height - 2):
            d2fdx = img[x+1, y] - 2*img[x,y] + img[x-1,y]
            d2fdy = img[x, y+1] - 2*img[x,y] + img[x,y-1]
            d2fdxdy = (img[x+1, y+1] + img[x-1, y-1] - img[x+1, y-1] - img[x-1, y+1])/4.0
            value = np.zeros((2, 2))
            value[0, 0] = d2fdx
            value[0, 1] = d2fdxdy
            value[1, 0] = d2fdxdy
            value[1, 1] = d2fdy
            hessianM[(x, y)] = value

    return hessianM


def getEigenValues(mat2x2):
    return np.linalg.eigvals(mat2x2)



def gradient(src):
    width = np.size(src, 0)
    height = np.size(src, 1)
    #gx = np.zeros((width, height))
    gx = createImg(height, width)
    gy = createImg(height, width)
    #gy = np.zeros((width, height))

    for x in range(1, width-2):
        for y in range(1, height-2):
            gx[x, y] = (float(src[x+1, y]) - float(src[x-1, y]))/2.0
            gy[x, y] = (float(src[x, y+1]) - float(src[x, y-1]))/2.0
            #if abs(gx[x,y]) <= 0.5 or abs(gy[x,y]) <= 0.5:
                #print(str(src[x+1, y]) + "-" + str(src[x-1, y]))
                #print(str(src[x, y+1]) + "-" + str(src[x, y-1]))

    return gx, gy


def magnitude(gx, gy):
    width = np.size(gx, 0)
    height = np.size(gx, 1)
    res = gx.copy()

    for x in range(1, width - 2):
        for y in range(1, height - 2):
            res[x, y] = math.sqrt(sqr(gx[x, y]) + sqr(gy[x, y]))

    return res


def sqr_magnitude(gx, gy):
    width = np.size(gx, 0)
    height = np.size(gx, 1)
    res = gx.copy()

    for x in range(1, width - 2):
        for y in range(1, height - 2):
            res[x, y] = sqr(gx[x, y]) + sqr(gy[x, y])

    return res


def almostEqual(v1, v2):
    if abs(v2 - v1) <= 0.00000000001:
        return True
    return False


def areLinearlyDependant(val1X, val2X, val1Y, val2Y):

    if val2X == 0 or val2Y == 0:
        if val1X == 0 and val1Y == 0 and val2X == 0 and val2Y == 0:
            return 0
    elif almostEqual(-val1X / val2X, -val1Y/val2Y):
        if includeNulLambda or not(almostEqual(-val1X / val2X, 0)):
            return 2

    if val1X == 0 or val1Y == 0:
        if val2X == 0 and val2Y == 0 and val1X == 0 and val1Y == 0:
            return 1
    elif almostEqual(-val2X / val1X, -val2Y / val1Y):
        if includeNulLambda or not(almostEqual(-val2X / val1X, 0)):
            return 2
    return -1


def getABCvalues(hessian):
    a = hessian[0, 0]
    b = hessian[0, 1]
    c = hessian[1, 1]

    return a, b, c

def isLambdaLinear(hessian, iGx, iGy):
    a, b, c = getABCvalues(hessian)
    if not(almostEqual(iGx, 0) and almostEqual(iGy, 0)):
        lambda1 = (a*iGx + b*iGy)/iGx
        lambda2 = (b*iGx + c*iGy)/iGy

        if almostEqual(lambda1, lambda2):
            return True

    elif (almostEqual(iGx, 0) and not(almostEqual(iGy, 0))) or (almostEqual(iGy, 0) and not(almostEqual(iGx, 0))):
        return almostEqual(b*(iGx+iGy), 0)

    else:
        return True

    return False


def getJacobi(f, g):
    width = np.size(f, 0)
    height = np.size(f, 1)
    res = list()
    resZeroG = list()
    resZeroSGM = list()

    fx, fy = gradient(f)
    gx, gy = gradient(g)
    for x in range(1, width - 2):
        for y in range(1, height - 2):
            linearity = areLinearlyDependant(fx[x, y], gx[x, y], fy[x, y], gy[x, y])
            if linearity == 2:
                res.append([x, y])

            if linearity == 1:
                resZeroG.append([x, y])

            if linearity == 0:
                resZeroSGM.append([x, y])

    return res, resZeroG, resZeroSGM


def pointsOnImg(src, pointlist, pointlistG, pointlistSGM):
    res = src.copy()
    for pt in pointlist:
        cv.circle(res, (pt[1], pt[0]), 2, (0,0,255),-1)

    if enableZeroG:
        for pt in pointlistG:
            cv.circle(res, (pt[1], pt[0]), 2, (0, 255, 0), -1)

        for pt in pointlistSGM:
            cv.circle(res, (pt[1], pt[0]), 2, (255, 0, 0), -1)

    return res


def computeCritical(listImg):
    for imgPath in listImg:
        print("Starting computation for "+imgPath)
        img = getImgGrayscale(path+imgPath)
        orig = cv.imread(path+imgPath)
        img = addNeumannBorder(img, 1)
        iGx, iGy = gradient(img)
        #iGm = magnitude(iGx, iGy)

        gfunc = sqr_magnitude(iGx, iGy)

        jacoList, jacoZeroG, jacoZeroSGM = getJacobi(img, gfunc)
        jMars = pointsOnImg(orig, jacoList, jacoZeroG, jacoZeroSGM)
        cv.imwrite(path+"cmp_criticals_"+imgPath, jMars)
        #cv.imwrite(path+"gray_"+imgPath, img)
        print(imgPath+" complete !")


def _main():
    global path
    global enableZeroG
    global includeNulLambda
    #path = "/home/sebastien/Documents/MUNI/IMGS/RV-Graph/"
    path = "./IMGS/"

    enableZeroG = False;
    includeNulLambda = True;

    imglist = ["generated.jpg"]
    #imglist.append("gray_binary_fingerprint.jpg")
    #imglist.append("gray_earth.jpg")
    #imglist.append("gray_brain.jpg")

    computeCritical(imglist)
    #cv.waitKey(0)
    #cv.destroyAllWindows()



_main()
