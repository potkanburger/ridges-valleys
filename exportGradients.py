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
        iGx, iGy = gradient(img)
        cv.imwrite(path + "x_gradient_" + imgPath, iGx)
        cv.imwrite(path + "y_gradient_" + imgPath, iGy)

        print(imgPath + " complete !")

def _main():
    global path
    # path = "/home/sebastien/Documents/MUNI/IMGS/RV-Graph/"
    path = "./IMGS/"

    imglist = ["generated.jpg"]
    # imglist.append("gray_binary_fingerprint.jpg")
    # imglist.append("gray_earth.jpg")
    # imglist.append("gray_brain.jpg")

    createGradients(imglist)


_main()