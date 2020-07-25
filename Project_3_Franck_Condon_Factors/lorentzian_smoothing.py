import numpy as np
from scipy.optimize import least_squares


def lorentzian(x, x0, a, gam):  # gam = 1/Γ
    return a * gam ** 2 / (gam ** 2 + (x - x0) ** 2)


def multi_lorentz(x, params):
    return sum([lorentzian(x, *params[i: i + 3]) for i in range(0, len(params), 3)])


def res_multi_lorentz(params, xData, yData):
    diff = [multi_lorentz(x, params) - y for x, y in zip(xData, yData)]
    return diff


def smooth(xData, yData, width=1, min_peaks=3, max_peaks=6, initial_guess_x=0, x0_step=1, ftol=1e-3, xtol=1e-3,
           gtol=1e-3):
    scale = max(yData)
    yData = yData / scale

    offset = min(yData)
    counter = 0

    residue = [y - offset for y in yData]
    optimization_result = None
    initial_params = []
    x0 = initial_guess_x
    while max(residue) > .1:
        counter += 1
        if counter > max_peaks:
            break
        max_res = np.argmax(residue)
        max_y = yData[max_res]
        initial_params += [x0, max_y - residue[max_res], width]
        optimization_result = least_squares(res_multi_lorentz, initial_params, args=(xData, yData), ftol=ftol,
                                            xtol=xtol,
                                            gtol=gtol, method='lm')
        residue = [y - multi_lorentz(x, optimization_result.x) for x, y in zip(xData, yData)]
        x0 += x0_step

    data = [offset + multi_lorentz(x, optimization_result.x) for x in xData]
    return data
