import cv2 as cv
import numpy as np
import math
from main import *


def createGradients(listImg):
    for imgPath in listImg:
        print("Starting computation for "+imgPath)
        img = getImgGrayscale(path+imgPath)
        orig = cv.imread(path+imgPath)
        img = addNeumannBorder(img, 1)

        width = np.size(img, 0)
        height = np.size(img, 1)

        iGx, iGy = gradient(img)
        cv.imwrite(path + "x_gradient_" + imgPath, iGx)
        cv.imwrite(path + "y_gradient_" + imgPath, iGy)
        cv.imwrite(path + "validator_" + imgPath, img)
        f = open(path + "x_gradient_" + imgPath +".txt", "w")

        for x in range(1, width - 2):
            for y in range(1, height - 2):
                f.write(str(iGx[x, y]))
            f.write("\n")
        f.close()
        print(imgPath + " complete !")

def _main():
    global path
    # path = "/home/sebastien/Documents/MUNI/IMGS/RV-Graph/"
    path = "./IMGS/"

    imglist = []
    # imglist.append("circular.jpg")
    imglist.append("smooth_lena.jpg")
    #imglist.append("ellipses_generated.jpg")
    # imglist.append("mars1.jpg")
    # imglist.append("gray_binary_fingerprint.jpg")
    # imglist.append("gray_earth.jpg")
    # imglist.append("gray_brain.jpg")

    createGradients(imglist)


_main()