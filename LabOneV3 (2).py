# Daniel Boe
# ME 553 - Laboratory Project One - Version Two
# January 18, 2019

from sympy import plot_implicit, symbols, Eq
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns


sns.set(style="white", rc={"axes.facecolor": (0, 0, 0, 0)})

# Define the system constants
h_0 = 5.75  #[m]
r_t = 3.25  #[m]
r_d = 0.015 #[m]
C_d = 0.6
g = 9.81  # [m/s^2]

# #This Section is used for plotting the Analytic Solution
# ##############################################################################

# # Decrlare the 'symbolic' variables for the implicit solver
# h, t = symbols('h t')

# # Integration Constant 
# phi = h_0**(3/2) *(3/10*h_0 - r_t)

# # Define the implicit equation
# eqn = Eq(h**(3/2)*(3/10*h - r_t) - 3/4*C_d *r_d**2 * (2*g)**(1/2)*t*3600,phi)
# plot = plot_implicit(eqn, x_var = (t, 0, 14), y_var = (h, 0, 6))

# # Access the plot in order to adjust Plots Appearance
# fig = plt.gcf()
# ax = plt.gca()

# # Create a grid for the plot
# #ax.grid()

# ax.set_xlabel(' ', horizontalalignment='center', x = 1.0, fontsize = 5)
# ax.set_ylabel(' ', verticalalignment = 'center', y = 1.0, fontsize = 5)

# ax.text(0.5, -0.1, 'Time [hours]', horizontalalignment='center', \
#           verticalalignment='center', transform=ax.transAxes, fontsize = 22)

# ax.text(-0.05, 0.5, 'Fluid Height [m]', horizontalalignment='center', \
#           verticalalignment='center', transform=ax.transAxes, fontsize = 22, rotation=90)

# ax.plot([0],[0], color = 'blue', label = 'Analytic\nSolution')
# ##############################################################################


# Set the number of points to discretize the fluid height domain
h_points = 1000

# Discretize the fluid height by creating an array that ranges from h_0 to 
# zero (use a small pertubation from zero to avoid division by zero).
hVals = np.linspace(h_0,0.000001,h_points)

# Calculate the Taylor Series constants for each h value
alpha = -r_d**2*C_d*np.sqrt(2*g)/(2*hVals**(1/2)*r_t - hVals**(3/2))
beta = r_d**2*C_d*np.sqrt(2*g)*(hVals**(-1/2)*r_t - 3/2*hVals**(1/2)) / \
    (2*hVals**(1/2)*r_t - hVals**(3/2))**2 

# Create an array to store all the time delta values
delta_t = np.zeros(len(hVals))

# Create a 'Shifted' h array which is essentially acting as the i'th element in 
# our discretized domain while the hVals acts as the (i-1) element.
h_i = np.append(hVals[1:],0)

# Solve for the time delta between each h value
delta_t = np.log((h_i - hVals)*beta/alpha + 1)/beta
delta_t2 = (h_i - hVals)/(alpha + beta*(h_i - hVals))

# Shift the delta_t array so that the first element in tVals is zero.
tVals = np.insert(delta_t[:-1],0,0)
tVals2 = np.insert(delta_t2[:-1],0,0)

# Perform an accumulating sum to create the time array
tVals = np.cumsum(tVals)
tVals2 = np.cumsum(tVals2)

# Create a Figure object
fig, ax = plt.subplots(1,1)
ax = [ax]
 
#ax[0].plot(tVals/3600, hVals, color = 'seagreen', label =\
#      'Linearized Model\nDiscretized Domain = ' + str(h_points) + ' points', linewidth = 6) 
ax[0].plot(tVals/3600, hVals, color = 'seagreen', label =\
      'Nonlinear Numerical Method\nDiscretized Domain = ' + str(h_points) + ' points', linewidth = 3)  
ax[0].plot(tVals2/3600, hVals, color = 'firebrick', label =\
      'Linear Numerical Method\nDiscretized Domain = ' + str(h_points) + ' points', linewidth = 3)  
    
# Set the axes titles    
ax[0].set_xlim(left = 0)
ax[0].set_xlim(right = 14)
ax[0].set_ylim(bottom = 0)
ax[0].set_ylim(top = 6)
ax[0].grid(linestyle = '-', linewidth = 2, alpha = 0.65)

# Creates overall legend for both plots
handles, labels = ax[0].get_legend_handles_labels()
ax[0].legend(handles, fontsize = 16, ncol =1, labels = labels, loc = 'upper right', labelspacing = 0.85, framealpha = 0.25)

#fig.suptitle('Linearization of $\dot{h}$ Centered About Various $h$ Values' , fontsize = 35)
fig.suptitle('Emptying Spherical Water Tank - ME 553: Laboratory Project One\nFluid Height vs. Time - Analytic and Linearized Solutions' \
             , fontsize = 22, weight = 'bold')

# Set the x label for x axis
#ax[0].set_xlabel('Time [hours]', fontsize=24)

# Set the y labels
#ax[0].set_ylabel('Fluid Height [meters]', fontsize = 24)


# Set the size of the axis labels
ax[0].tick_params(axis='y', which='major', labelsize=18)
plt.tick_params(axis = 'x', which = 'major', labelsize = 18,  direction = 'out')

plt.show()