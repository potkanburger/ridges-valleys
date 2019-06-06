import cv2 as cv
import numpy as np
import math
from main import *

def getHeightPointsList(img):
    iGx, iGy = gradient(img)
    resHeight = list()

    width = np.size(img, 0)
    height = np.size(img, 1)

    hessianMatrix = hessian(img)

    for x in range(1, width - 2):
        for y in range(1, height - 2):
            if not(almostEqual(iGx[x, y], 0) and almostEqual(iGy[x,y], 0)):
                ev = getEigenValues(hessianMatrix[(x, y)])
                if not(almostEqual(ev[0], 0) or almostEqual(ev[1], 0) or almostEqual(ev[0], ev[1])):
                    print(ev)
                    if isLambdaLinear(hessianMatrix[(x, y)], iGx[x, y], iGy[x,y]):
                        resHeight.append((x, y))
    return resHeight

def computeHeightPoints(listImg):
    for imgPath in listImg:
        print("Starting height point computation for "+imgPath)
        img = getImgGrayscale(path+imgPath)
        orig = cv.imread(path+imgPath)
        hList = getHeightPointsList(img)
        hpImg = pointsOnImg(orig, hList, list(), list())
        cv.imwrite(path+"hpts_"+imgPath, hpImg)
        print(imgPath+" complete !")


def _main():
    global path
    path = "./IMGS/"
    imglist = []
    imglist.append("ellipses_generated_rotated90.jpg")
    #imglist.append("ellipses_generated.jpg")
    #imglist.append("generated.jpg")
    #imglist.append("gray_brain.jpg")
    # imglist.append("gray_binary_fingerprint.jpg")
    # imglist.append("gray_earth.jpg")
    # imglist.append("gray_brain.jpg")
    computeHeightPoints(imglist)


_main()