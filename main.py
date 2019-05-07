import cv2 as cv
import numpy as np
import math
import matplotlib.pyplot as plot

def getImgGrayscale(url):
    return cv.imread(url, 0)

def sqr(val):
    return pow(val, 2)

def gradient(src):
    width = np.size(src, 0)
    height = np.size(src, 1)
    #gx = np.zeros((width, height))
    gx = src.copy()
    gy = src.copy()
    #gy = np.zeros((width, height))

    for x in range(1, width-2):
        for y in range(1, height-2):
            gx[x, y] = (int(src[x+1, y]) - int(src[x-1, y]))/2.0
            gy[x, y] = (int(src[x, y+1]) - int(src[x, y-1]))/2.0

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
    if(abs(v2-v1) <= 0.00000000001):
        return True
    return False

def areLinearlyDependant(val1X, val2X, val1Y, val2Y):

    if val2X == 0 or val2Y == 0:
        if val1X == 0 and val1Y == 0 and val2X == 0 and val2Y == 0:
            return 1
    elif almostEqual(-val1X / val2X, -val1Y/val2Y):
        return 2

    if val1X == 0 or val1Y == 0:
        if val2X == 0 and val2Y == 0 and val1X == 0 and val1Y == 0:
            return 1
    elif almostEqual(-val2X / val1X, -val2Y/ val1Y):
        return 2

    return 0


def getJacobi(f, g):
    width = np.size(f, 0)
    height = np.size(f, 1)
    res = list()
    resZeroG = list()
    fx, fy = gradient(f)
    gx, gy = gradient(g)
    for x in range(1, width - 2):
        for y in range(1, height - 2):
            linearity = areLinearlyDependant(fx[x, y],gx[x, y], fy[x, y], gy[x, y])
            if linearity == 2:
                res.append([x, y])

            if linearity == 1:
                resZeroG.append([x, y])

    return res, resZeroG

def pointsOnImg(src, pointlist, pointlist2):
    res = src.copy()
    for pt in pointlist:
        cv.circle(res, (pt[1], pt[0]), 2, (0,0,255),-1)

    for pt in pointlist2:
        cv.circle(res, (pt[1], pt[0]), 2, (0, 255, 0), -1)

    return res

def computeCritical(listImg):
    for imgPath in listImg:
        print("Starting computation for "+imgPath)
        img = getImgGrayscale(path+imgPath)
        orig = cv.imread(path+imgPath)
        iGx, iGy = gradient(img)
        iGm = magnitude(iGx, iGy)

        gfunc = sqr_magnitude(iGx, iGy)

        jacoList, jacoZeroG = getJacobi(img, gfunc)
        jMars = pointsOnImg(orig, jacoList, jacoZeroG)
        cv.imwrite(path+"cmp_criticals_"+imgPath, jMars)
        #cv.imwrite(path+"gray_"+imgPath, img)
        print(imgPath+" complete !")

def main():
    global path
    path = "/home/sebastien/Documents/MUNI/IMGS/RV-Graph/"

    imglist = ["gray_mars1.jpg"]
    imglist.append("gray_binary_fingerprint.jpg")
    imglist.append("gray_earth.jpg")
    imglist.append("gray_brain.jpg")

    computeCritical(imglist)
    #cv.waitKey(0)
    #cv.destroyAllWindows()



main()
