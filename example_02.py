'''
Attempt to implement a non-linear pendulum with perturbations using:

    alpha = - g / l * sin (theta) + (input_force + damping_force) / m / l

This scripts generates and plots the behavior of `theta` with respect to time and a rudimental
animation of the pendulum.

The initial condition can be assigned by setting `forced_torque_vec` as a vector of torque applied
during the first N time steps

'''

import time
import numpy as np
import matplotlib.pyplot as plt

g = 9.81  # m / s^2
dt = 0.001  # s

######################################## Parameters ###############################################
animation = True
DURATION = 10  # s
l = 0.5  # m
m = 0.5  # kg
# forced_torque_vec = [70 for _ in range(0, 60)]
forced_torque_vec = [2 * np.sin(np.sqrt(g / l) * i * dt) for i in range(0, 4000)]
damping_k = 0.8  # kg / s^2
#################################### End of parameters ############################################

NUM_OF_SAMPLES = round(DURATION / dt)


def myPlotter(xVals, yVals, xLabel, yLabel, title):
    plt.figure(figsize=(8, 6))
    plt.plot(xVals, yVals, color='blue', linewidth=4)
    plt.title(title)
    plt.xlabel(xLabel)
    plt.ylabel(yLabel)
    plt.grid(visible=True)


def theta_t_plus_dt(dt, theta_t, theta_t_minus_dt, forced_torque=0, damping_k=0):
    alpha_f = damping_k * (theta_t_minus_dt - theta_t) / dt / m
    alpha_in = forced_torque / l / l / m
    alpha_g = -g / l * np.sin(theta_t)
    alpha = alpha_g + alpha_in + alpha_f
    return alpha * (dt * dt) - theta_t_minus_dt + 2 * theta_t


t = np.linspace(0, NUM_OF_SAMPLES * dt, NUM_OF_SAMPLES + 1)
theta_curr = 0
theta_prev = 0
theta_next = 0
theta_vec = [0 for _ in range(0, NUM_OF_SAMPLES + 1)]
omega_vec = [0 for _ in range(0, NUM_OF_SAMPLES + 1)]
torque_vec = [0 for _ in range(0, NUM_OF_SAMPLES + 1)]
len_of_forced_torque_vec = min(NUM_OF_SAMPLES + 1, len(forced_torque_vec))

# If there is any user-defined torque, this will be used as a system input
for i in range(0, len_of_forced_torque_vec):
    torque_vec[i] = forced_torque_vec[i]


start_time = time.time()
for i in range(0, NUM_OF_SAMPLES + 1):
    theta_next = theta_t_plus_dt(dt, theta_curr, theta_prev, torque_vec[i], damping_k)
    omega_vec[i] = (theta_next - theta_curr) / dt
    theta_prev = theta_curr
    theta_curr = theta_next
    theta_vec[i] = theta_curr
print("Computation time =", time.time() - start_time)

myPlotter(t, theta_vec, 't [s]', 'theta [rad]', 'angle vs time')
myPlotter(t, torque_vec, 't [s]', 'applied torque [Nm]', 'Applied torque')

# Energy / state space
myPlotter(theta_vec, omega_vec, 'theta [rad]', 'omega [rad/s]', 'State Space')
KE = [(omega_vec[i] * l)**2 * m / 2 for i in range(0, NUM_OF_SAMPLES + 1)]
U = [m * g * l * (1 - np.cos(theta_vec[i])) for i in range(0, NUM_OF_SAMPLES + 1)]
E = [KE[i] + U[i] for i in range(0, NUM_OF_SAMPLES + 1)]
myPlotter(KE, U, 'kinetic [J]', 'potential [J]', 'Energy')
myPlotter(t, E, 'time [s]', 'total energy [J]', 'Energy vs time')

if (animation):
    fig, ax = plt.subplots()
    ax.set_xlim(-l * 1.1, l * 1.1)
    ax.set_ylim(-l * 1.1, l * 1.1)
    ax.set_aspect('equal', adjustable='box')
    position_plot, = ax.plot(0, 0, marker='*')
    refresh_rate = 50  # Hz
    undersampling_rate = round(1 / dt / refresh_rate)
    for i in range(0, NUM_OF_SAMPLES + 1, undersampling_rate):
        position_plot.set_data([0, l * np.sin(theta_vec[i])], [0, - l * np.cos(theta_vec[i])])
        plt.pause(1 / refresh_rate)

plt.show()
