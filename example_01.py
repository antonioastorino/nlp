'''
Attempt to implement a non-linear pendulum using:

    alpha = - g / l * sin (theta)

This scripts generates and plots the behavior of `theta` with respect to time and a rudimental
animation of the pendulum.

> NOTE: I cannot find a better implementation online so I made this up and I am not sure this is
        correct. It seems that the energy is conserved. To see that, set
        * `animation = False`
        * `DURATION = 1000`

> NOTE: python is slow at rendering and your frame rate is anyway the bottleneck. For this reason,
    the `dt` used for generating data is not the same used to update the animation
'''

import numpy as np
import matplotlib.pyplot as plt

######################################## Parameters ###############################################
animation = True
INITIAL_THETA = np.pi * 3 / 4
DURATION = 6  # s
l = 2  # m
m = 10  # kg
#################################### End of parameters ############################################

dt = 0.001  # s
g = 9.81  # m / s^2
NUM_OF_SAMPLES = round(DURATION / dt)


def myPlotter(xVals, yVals, xLabel, yLabel, title):
    plt.figure(figsize=(8, 6))
    plt.plot(xVals, yVals, color='blue', linewidth=4)
    plt.title(title)
    plt.xlabel(xLabel)
    plt.ylabel(yLabel)
    plt.grid(visible=True)


def theta_t_plus_dt(dt, theta_t, theta_t_minus_dt):
    return -g / l * np.sin(theta_t) * (dt * dt) - theta_t_minus_dt + 2 * theta_t


t = np.linspace(0, NUM_OF_SAMPLES * dt, NUM_OF_SAMPLES - 1)
theta_curr = INITIAL_THETA
theta_prev = INITIAL_THETA
theta_next = 0
theta_vec = [0 for _ in range(0, NUM_OF_SAMPLES - 1)]
omega_vec = [0 for _ in range(0, NUM_OF_SAMPLES - 1)]
for i in range(0, NUM_OF_SAMPLES - 1):
    theta_next = theta_t_plus_dt(dt, theta_curr, theta_prev)
    omega_vec[i] = (theta_next - theta_curr) / dt
    theta_prev = theta_curr
    theta_curr = theta_next
    theta_vec[i] = theta_curr

myPlotter(t, theta_vec, 't [s]', 'theta [rad]', 'angle vs time')

# Energy / state space
myPlotter(theta_vec, omega_vec, 'theta', 'omega', 'State Space')
KE = [(omega_vec[i] * l)**2 * m / 2 for i in range(0, NUM_OF_SAMPLES - 1)]
U = [m * g * l * (1 - np.cos(theta_vec[i])) for i in range(0, NUM_OF_SAMPLES - 1)]
E = [KE[i] + U[i] for i in range(0, NUM_OF_SAMPLES - 1)]
myPlotter(KE, U, 'kinetic [J]', 'potential [J]', 'Energy')
myPlotter(t, E, 'time [s]', 'total energy [J]', 'Energy vs time')
print("------------------------------------")
print("Energy deviation from average: {:3.2f}%".format(
    100 * (max(E) - min(E)) / (max(E) + min(E) * 2)))
print("------------------------------------")

if (animation):
    fig, ax = plt.subplots()
    ax.set_xlim(-l * 1.1, l * 1.1)
    ax.set_ylim(-l * 1.1, l * 1.1)
    ax.set_aspect('equal', adjustable='box')
    position_plot, = ax.plot(0, 0, marker='*')
    refresh_rate = 50  # Hz
    undersampling_rate = round(1 / dt / refresh_rate)
    for i in range(0, NUM_OF_SAMPLES - 1, undersampling_rate):
        position_plot.set_data([0, l * np.sin(theta_vec[i])], [0, - l * np.cos(theta_vec[i])])
        plt.pause(1 / refresh_rate)

plt.show()
