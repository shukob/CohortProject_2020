import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import leastsq


def lorentzian(x, x0, a, gam):
    return a * gam ** 2 / (gam ** 2 + (x - x0) ** 2)


def multi_lorentz(x, params):
    off = params[0]
    paramsRest = params[1:]
    assert not (len(paramsRest) % 3)
    return off + sum([lorentzian(x, *paramsRest[i: i + 3]) for i in range(0, len(paramsRest), 3)])


def res_multi_lorentz(params, xData, yData):
    diff = [multi_lorentz(x, params) - y for x, y in zip(xData, yData)]
    return diff


def smooth(xData, yData, width=1):
    yData = yData / max(yData)

    width = 1

    yDataLoc = yData
    startValues = [max(yData)]
    counter = 0

    while max(yDataLoc) - min(yDataLoc) > .1:
        counter += 1
        if counter > 20:  ### max 20 peak...emergency break to avoid infinite loop
            break
        minP = np.argmin(yDataLoc)
        minY = yData[minP]
        x0 = xData[minP]
        startValues += [x0, minY - max(yDataLoc), width]
        popt, ier = leastsq(res_multi_lorentz, startValues, args=(xData, yData))
        yDataLoc = [y - multi_lorentz(x, popt) for x, y in zip(xData, yData)]

    return yDataLoc

