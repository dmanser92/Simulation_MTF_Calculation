# -*- coding: utf-8 -*-
"""
Created on Wed Sep 27 15:29:01 2017

@author: axdma
"""
# importieren von Bibliotheken
import matplotlib.pyplot as plt
import numpy as np
import scipy.optimize


# definiern von Funktionen mit dem lambda-Operator
def sin(x):
    return np.sin(np.deg2rad(x))


# fit a sine to a matrix
def fit_sin(tt, yy):
    # fit sin to the input time sequence, and return fitting parameters
    # "amp", "omega", "phase", "offset", "freq", "period" and "fitfunc"'''
    tt = np.array(tt)
    ff = np.fft.fftfreq(len(tt), (tt[1] - tt[0]))  # assume uniform spacing
    yy_stored = yy
    guess_freq = 0
    guess_amp = 0
    guess_offset = 0

    for i in range(0, len(yy_stored)):
        yy = np.array(yy_stored[i])
        Fyy = abs(np.fft.fft(yy))
        guess_freq += abs(ff[np.argmax(Fyy[1:])+1])   # excluding the zero frequency "peak", which is related to offset
        guess_amp += np.std(yy) * 2.**0.5
        guess_offset += np.mean(yy)

    guess = np.array([guess_amp, 2. * np.pi * guess_freq, 0., guess_offset]) / len(yy_stored)

    def sinfunc(t, A, w, p, c): return A * np.sin(w*t + p) + c
    popt, pcov = scipy.optimize.curve_fit(sinfunc, tt, yy, p0=guess)
    A, w, p, c = popt
    f = w/(2.*np.pi)
    fitfunc = lambda t: A * np.sin(w*t + p) + c
    return {"amp": A, "omega": w, "phase": p, "offset": c, "freq": f, "period": 1./f, "fitfunc": fitfunc, "maxcov": np.max(pcov), "rawres": (guess,popt,pcov)}


N = 50
t = np.linspace(0, 10, N)
t2 = np.linspace(0, 10, 10*N)
F = 4*sin(2*np.pi*8*t+0.4) + 5
intensityValues = np.zeros((N, N))

for i in range(0, N):
    for j in range(0, N):
        intensityValues[i][j] = F[j] + 2 * np.random.sample() - 1   # continuous uniform distribution
#        intensityValues[i][j] = F[j] + np.random.normal(0, 0.3)    # normal distribution (gaussian)


res = fit_sin(t, intensityValues)
print("Amplitude=%(amp)s, Angular freq.=%(omega)s, phase=%(phase)s, offset=%(offset)s, Max. Cov.=%(maxcov)s" % res)

plt.plot(t, F, "-k", label="y", linewidth=2)
plt.plot(t, intensityValues.T, "oy", label="y with noise")
plt.plot(t2, res["fitfunc"](t2), "r-", label="y fit curve", linewidth=2)
plt.legend(loc="best")
plt.xlabel('Pixelachse x [mm]', fontsize=12)
plt.ylabel('Intensitätswerte der Pixel [W/m2]', fontsize=12)
plt.title('Simulation MTF Messpunkte', fontsize=20)
plt.show()
