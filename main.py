import cv2 as cv
import numpy as np
import math
#import matplotlib.pxplot as plot

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

def scale0to250(img):
    with np.nditer(img, op_flags=['readwrite']) as it:
        for x in it:
            x[...] = x+125.0

def hessian(img):
    width = np.size(img, 0)
    height = np.size(img, 1)
    hessianM = {}
    for y in range(1, width - 2):
        for x in range(1, height - 2):
            d2fdy = img[y+1, x] - 2*img[y,x] + img[y-1,x]
            d2fdx = img[y, x+1] - 2*img[y,x] + img[y,x-1]
            d2fdydx = (img[y+1, x+1] + img[y-1, x-1] - img[y+1, x-1] - img[y-1, x+1])/4.0
            value = np.zeros((2, 2))
            value[0, 0] = d2fdx
            value[0, 1] = d2fdydx
            value[1, 0] = d2fdydx
            value[1, 1] = d2fdy
            hessianM[(y, x)] = value

    return hessianM


def getEigenValues(mat2x2):
    return np.linalg.eigvals(mat2x2)



def gradient(src):
    width = np.size(src, 0)
    height = np.size(src, 1)
    #gy = np.zeros((width, height))
    gy = createImg(height, width)
    gx = createImg(height, width)
    #gx = np.zeros((width, height))

    for y in range(1, width-2):
        for x in range(1, height-2):
            gx[y, x] = (float(src[y, x+1]) - float(src[y, x-1]))/2.0
            gy[y, x] = (float(src[y + 1, x]) - float(src[y - 1, x])) / 2.0
            #if abs(gy[x,y]) <= 0.5 or abs(gx[x,y]) <= 0.5:
                #print(str(src[x+1, y]) + "-" + str(src[x-1, y]))
                #print(str(src[x, y+1]) + "-" + str(src[x, y-1]))

    return gx, gy


def magnitude(gy, gx):
    width = np.size(gy, 0)
    height = np.size(gy, 1)
    res = gy.copy()

    for y in range(1, width - 2):
        for x in range(1, height - 2):
            res[y, x] = math.sqrt(sqr(gy[y, x]) + sqr(gx[y, x]))

    return res


def sqr_magnitude(gy, gx):
    width = np.size(gy, 0)
    height = np.size(gy, 1)
    res = gy.copy()

    for y in range(1, width - 2):
        for x in range(1, height - 2):
            res[y, x] = sqr(gy[y, x]) + sqr(gx[y, x])

    return res


def almostEqual(v1, v2):
    if abs(v2 - v1) <= 0.00000000001:
        return True
    return False


def areLinearlyDependant(val1y, val2y, val1x, val2x):

    if val2y == 0 or val2x == 0:
        if val1y == 0 and val1x == 0 and val2y == 0 and val2x == 0:
            return 0
    elif almostEqual(-val1y / val2y, -val1x/val2x):
        if includeNulLambda or not(almostEqual(-val1y / val2y, 0)):
            return 2

    if val1y == 0 or val1x == 0:
        if val2y == 0 and val2x == 0 and val1y == 0 and val1x == 0:
            return 1
    elif almostEqual(-val2y / val1y, -val2x / val1x):
        if includeNulLambda or not(almostEqual(-val2y / val1y, 0)):
            return 2
    return -1


def getABCvalues(hessian):
    a = hessian[0, 0]
    b = hessian[0, 1]
    c = hessian[1, 1]

    return a, b, c

def isLambdaLinear(hessian, iGy, iGx):
    a, b, c = getABCvalues(hessian)
    if not(almostEqual(iGy, 0)) and not(almostEqual(iGx, 0)):
        lambda1 = (a*iGy + b*iGx)/iGy
        lambda2 = (b*iGy + c*iGx)/iGx

        if almostEqual(lambda1, lambda2):
            return True

    elif (almostEqual(iGy, 0) and not(almostEqual(iGx, 0))) or (almostEqual(iGx, 0) and not(almostEqual(iGy, 0))):
        return almostEqual(b*(iGy+iGx), 0)

    else:
        return True

    return False


def getJacobi(f, g):
    width = np.size(f, 0)
    height = np.size(f, 1)
    res = list()
    resZeroG = list()
    resZeroSGM = list()

    fy, fx = gradient(f)
    gy, gx = gradient(g)
    for y in range(1, width - 2):
        for x in range(1, height - 2):
            linearity = areLinearlyDependant(fy[y, x], gy[y, x], fx[y, x], gx[y, x])
            if linearity == 2:
                res.append([y, x])

            if linearity == 1:
                resZeroG.append([y, x])

            if linearity == 0:
                resZeroSGM.append([y, x])

    return res, resZeroG, resZeroSGM


def pointsOnImg(src, pointlist, pointlistG, pointlistSGM):
    res = src.copy()
    for pt in pointlist:
        cv.circle(res, (pt[0], pt[1]), 2, (0,0,255),-1)

    if enableZeroG:
        for pt in pointlistG:
            cv.circle(res, (pt[0], pt[1]), 2, (0, 255, 0), -1)

        for pt in pointlistSGM:
            cv.circle(res, (pt[0], pt[1]), 2, (255, 0, 0), -1)

    return res


def computeCritical(listImg):
    for imgPath in listImg:
        print("Starting computation for "+imgPath)
        img = getImgGrayscale(path+imgPath)
        orig = cv.imread(path+imgPath)
        img = addNeumannBorder(img, 1)
        iGy, iGx = gradient(img)
        #iGm = magnitude(iGy, iGx)

        gfunc = sqr_magnitude(iGy, iGx)

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
