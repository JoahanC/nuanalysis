import numpy as np 
import matplotlib.pyplot as plt
from lc_gen import *
from scipy.optimize import curve_fit


def func(x, a, b, c):
    return a * np.exp(-b * x) + c


def lightcurve_fit(x_data, y_data, char_time):
    """
    Applies a power law fit as a function of time elapsed from the characteristic time of the 
    peak of a transient pulse.

    Arguments
    ---------
    dataset : `array-like`
        A two dimensional array containing lightcurve data
    peak : `float`
        A value specifying the peak of the lightcurve
    char_time : `float`
        A value specifying the minimum reasonable time for a power law fit.
    """

    x_data = np.array(x_data)
    y_data = np.array(y_data)
    t_elapse = np.max(x_data) - char_time
    if t_elapse < 5:
        raise ValueError("Remaining curve fitting space is not large enough!")
    
    val_x_data = x_data
    val_y_data = y_data
    end_cuts = np.arange(np.min(val_x_data) + char_time, np.max(val_x_data), 0.1)
    fit_intervals = []
    for time in end_cuts:
        trim_x = val_x_data[val_x_data < time]
        trim_y = val_y_data[:len(trim_x)]
        fit_intervals.append([trim_x, trim_y])

    fit_vals = []
    for interval in fit_intervals:
        popt, popc = curve_fit(func, interval[0], interval[1], bounds=(0, [np.max(interval[1]) + 4, 5, interval[1][-1] + 4]))
        fit_vals.append([popt, popc])

    return end_cuts, fit_vals


time_domain = np.linspace(0, 1000, 10000)
data, data_err = gen_multicurve(time_domain, 5, 3, 100)

fig, ax = plt.subplots()
ax.errorbar(time_domain, data, data_err, ls='none', fmt="b+", ecolor="red")
ax.set_xlabel("Time")
ax.set_ylabel("Counts/seconds")
ax.set_title("Simulated NuSTAR Data")
plt.show()

max_idx = np.argmax(data)
fig, ax = plt.subplots()
if 10000 - max_idx >= 205:
    fit_time = time_domain[max_idx: max_idx + 200] - time_domain[max_idx]
    fit_data = data[max_idx: max_idx + 200]
    fit_err = data_err[max_idx: max_idx + 200]
    ax.errorbar(fit_time, fit_data, fit_err, ls='none', fmt="b+", ecolor="red")
    time_cuts, fit_vals = lightcurve_fit(fit_time, fit_data, 7)
if 10000 - max_idx < 205:
    fit_time = time_domain[max_idx:] - time_domain[max_idx]
    fit_data = data[max_idx:]
    fit_err = data_err[max_idx:]
    ax.errorbar(fit_time, fit_data, fit_err, ls='none', fmt="b+", ecolor="red")
    time_cuts, fit_vals = lightcurve_fit(fit_time, fit_data, 7)
ax.set_xlabel("Time")
ax.set_ylabel("Counts/seconds")
ax.set_title("Simulated NuSTAR Data")
plt.show()


fig, ax = plt.subplots()
ax.errorbar(fit_time, fit_data, fit_err, ls='none', fmt="b+", ecolor="red")
for val in fit_vals:
    amp = val[0][0]
    decay = val[0][1]
    offset = val[0][2]
    ax.plot(fit_time, func(fit_time, amp, decay, offset), color="black")

ax.set_xlabel("Time")
ax.set_ylabel("Counts/sec")
plt.show()

fig, ax = plt.subplots() 
gamma_vals = []
for val in fit_vals:
    gamma_vals.append(val[0][1])
ax.scatter(time_cuts, gamma_vals, s=2, color="black")
ax.set_xlabel("Time")
ax.set_ylabel(r"$\Gamma$")
ax.set_ylim(min(0.5, np.min(gamma_vals)), max(2.5, np.max(gamma_vals)))
plt.show()