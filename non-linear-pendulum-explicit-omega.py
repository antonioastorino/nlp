'''
Attempt to implement a non-linear pendulum using:

    alpha = - g / l * sin (theta)

This scripts generates and plots the behavior of `theta` with respect to time and a rudimental
animation of the pendulum.

> NOTE: I cannot find a better implementation online so I made this up and I am not sure this is
        correct. It seems that the energy is conserved. To see that, set
        * `animation = False`
        * `DURATION = 1000`

The problem is setting the initial conditions. Basically, you can define them in terms of initial
displacement (theta(t = -dt) and theta(t = 0)). This is causes problems if you change the
time interval `dt`, as the same displacement would occur within a different time slot. So, don't
change `dt` or try to also change the initial delta `theta` accordingly.

> NOTE: python is slow at rendering and your frame rate is anyway the bottleneck. For this reason,
    the `dt` used for generating data is not the same used to update the animation
> NOTE2: the simulation works also if the pendulum starts spinning around but wrapping `theta`
        when it exceeds 2pi makes the calculation of the angular speed `omega` more complicated.
        So, the energy-related plots will look weird. Test by using
        * `INITIAL_DELTA_THETA = np.pi / 530`
'''

import numpy as np
import matplotlib.pyplot as plt

######################################## Parameters ###############################################
animation = True
INITIAL_ENERGY = 370
DURATION = 10  # s
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


t = np.linspace(0, NUM_OF_SAMPLES * dt, NUM_OF_SAMPLES - 1)
theta_vec = [0 for _ in range(0, NUM_OF_SAMPLES - 1)]
omega_vec = [0 for _ in range(0, NUM_OF_SAMPLES - 1)]
alpha_vec = [0 for _ in range(0, NUM_OF_SAMPLES - 1)]
KE = [0 for _ in range(0, NUM_OF_SAMPLES - 1)]
U = [0 for _ in range(0, NUM_OF_SAMPLES - 1)]
E = [0 for _ in range(0, NUM_OF_SAMPLES - 1)]
theta_curr = 0
omega_sign = 1

for i in range(0, NUM_OF_SAMPLES - 1):
    alpha_vec[i] = -g / l * np.sin(theta_curr)
    omega_vec[i] = omega_sign * np.sqrt(2 * INITIAL_ENERGY / m / l / l - 2 *
                                        g / l * (1 - np.cos(theta_curr)))
    theta_vec[i] = theta_curr + omega_vec[i] * dt + 1 / 2 * alpha_vec[i] * dt * dt
    KE[i] = (omega_vec[i] * l)**2 * m / 2
    U[i] = m * g * l * (1 - np.cos(theta_vec[i]))
    E[i] = KE[i] + U[i]
    theta_curr = theta_vec[i]
    if i == 0:
        continue
    if np.abs(KE[i]) < 0.0001 and KE[i] < KE[i - 1]:
        omega_sign = -omega_sign
        print(KE[i])

myPlotter(t, theta_vec, 't [s]', 'theta [rad]', 'angle vs time')
myPlotter(t, omega_vec, 't [s]', 'omega [rad/s]', 'velocity vs time')
myPlotter(t, alpha_vec, 't [s]', 'alpha [rad/s/s]', 'acceleration vs time')

# Energy / state space
myPlotter(theta_vec, omega_vec, 'theta', 'omega', 'State Space')
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
