import cv2 as cv
import numpy as np
import math
from main import *


def _main():

    sizeDef = 300
    img = createImg(sizeDef, sizeDef)

    maxPoints = [(150, 150), (250, 250)]

    for x in range(sizeDef):
        for y in range(sizeDef):
            distance = sizeDef
            for pt in maxPoints:
                tmp_dst = math.sqrt(sqr(x - pt[0]) + sqr((y - pt[1])/2))
                if tmp_dst < distance:
                    distance = tmp_dst

            value = float(255 - distance)
            img[y, x] = max(float(0.0), value)

    cv.imwrite("./IMGS/ellipses_generated.jpg", img)


_main()
