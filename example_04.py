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
import random

g = 9.81  # Gravitational acceleration [m / s^2]
dt = 0.001  # Time step [s]

######################################## Parameters ###############################################
animation = True
DURATION = 10  # Simulation duration [s]
l = 1  # Pendulum length [m]
m = 2  # [kg]
damping_k = 0.8  # Damping coefficient [kg / s^2]
theta_d = np.pi / 4 * 4  # Desired angle [rad]
TORQUE_TO_CURR = 4  # Torque to current ratio [Nm/A]
MAX_CURRENT = 15  # Max controlling current [A]
K_p = 200
K_i = 10
K_d = 20
CURRENT_TF_POLE = 10
perturbation = [0 for _ in range(2000)] + [-200 for _ in range(10)]
SIGMA_MEAS = 0.01 * np.pi
SIGMA_CTRL = MAX_CURRENT * 0.001
#################################### End of parameters ############################################

NUM_OF_SAMPLES = round(DURATION / dt)


def myPlotter(xVals, yVals, xLabel, yLabel, title):
    plt.figure(figsize=(8, 6))
    plt.plot(xVals, yVals, color='blue', linewidth=2)
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


def i_t_plus_dt(dt, i_t, i_d):
    i_prime = i_t + (CURRENT_TF_POLE * (i_d - i_t)) * dt
    return MAX_CURRENT * np.tanh(i_prime / MAX_CURRENT)


t = np.linspace(0, NUM_OF_SAMPLES * dt, NUM_OF_SAMPLES + 1)
theta_curr = 0
theta_prev = 0
theta_next = 0
current_curr = 0
current_next = 0
theta_vec = [0 for _ in range(0, NUM_OF_SAMPLES + 1)]
omega_vec = [0 for _ in range(0, NUM_OF_SAMPLES + 1)]
torque_vec = [0 for _ in range(0, NUM_OF_SAMPLES + 1)]
current_vec = [0 for _ in range(0, NUM_OF_SAMPLES + 1)]
theta_err = [0 for _ in range(0, NUM_OF_SAMPLES + 1)]
perturbation_vec = [0 for _ in range(0, NUM_OF_SAMPLES + 1)]
perturbation_vec[0:min(NUM_OF_SAMPLES + 1, len(perturbation))] = perturbation[0:]
measurement_noise = [random.gauss(0, SIGMA_MEAS) for _ in range(0, NUM_OF_SAMPLES + 1)]
control_noise = [random.gauss(0, SIGMA_CTRL) for _ in range(0, NUM_OF_SAMPLES + 1)]
err_integral = 0
err_derivative = 0

start_time = time.time()
for i in range(0, NUM_OF_SAMPLES + 1):
    theta_err[i] = (theta_d - theta_curr + measurement_noise[i])
    err_integral = err_integral + theta_err[i]
    if i > 0:
        err_derivative = theta_err[i] - theta_err[i - 1]
    forced_current = K_p * theta_err[i] + K_i * (err_integral * dt) + K_d * (err_derivative / dt)
    current_next = i_t_plus_dt(dt, current_curr, forced_current)
    ctrl_torque = TORQUE_TO_CURR * (current_next + control_noise[i])
    theta_next = theta_t_plus_dt(
        dt,
        theta_curr,
        theta_prev,
        ctrl_torque +
        perturbation_vec[i],
        damping_k)
    omega_vec[i] = (theta_next - theta_curr) / dt
    theta_prev = theta_curr
    theta_curr = theta_next
    current_curr = current_next
    theta_vec[i] = theta_curr
    torque_vec[i] = ctrl_torque
    current_vec[i] = current_curr

print("Computation time =", time.time() - start_time)
print("Final error % =", (theta_d - theta_vec[-1]) / theta_d * 100)

theta_plot_vec = [theta_vec[i] / np.pi for i in range(0, NUM_OF_SAMPLES + 1)]
err_plot_vec = [theta_err[i] / theta_d * 100 for i in range(0, NUM_OF_SAMPLES + 1)]
myPlotter(t, theta_plot_vec, 't [s]', r"$\frac{\theta}{\pi}$ [rad]", 'angle vs time')
myPlotter(t, err_plot_vec, 't [s]', r"$(\theta_d - \theta)/\theta_d * 100$[%]", 'error vs time')
myPlotter(t, torque_vec, 't [s]', r'$\tau$ [Nm]', 'Controlling torque')
myPlotter(t, current_vec, 't [s]', r'i [A]', 'Controlling current')
myPlotter(theta_vec, omega_vec, 'theta [rad]', 'omega [rad/s]', 'State Space')

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
