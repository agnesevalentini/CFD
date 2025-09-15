"""
Compute boundary-layer thickness (delta_0.xx) for a set of CSVs and plot them.

Usage:
- Put this script in a folder with subfolder 'csvs' or change `csv_folder`.
"""

import numpy as np
import csv
import matplotlib.pyplot as plt

# ----- USER SETTINGS -----
x = [0.03, 0.06, 0.1, 0.2, 0.3, 0.4] # vector of the points on the foil chosen to extract the delta
chord = 0.48 # chord of the foil
t = 0.12 # thickness of the foil
Re = 4000
csv_folder = 'airfoil_{}_noSlip_refined'   # folder containing the csv files, {} -> Re
csv_name = 'u_airfoil_refined_x{}_{}.csv' # the {} will be substituted by x and Re
Bl_thickness = 0.99 # thickness of boundary layer, it should be between 0 and 1
slip = True # if true, the BL is computed using the U_e extracted from the slip case
slip_csv_folder = 'airfoil_{}_slip_refined' # folder of the slip csv file, ignored if slip = False
slip_csv_name = 'u_airfoil_refined_x{}_{}_slip.csv' # name of slip csv file, ignored if slip = False
# -------------------------

# We need this function to compute the normal and the tangent to the foil
def foil_derivative(t:float, chord:float, x):
    xp = x/chord
    foil_der = 5*t*(0.2969*np.power(xp,-1/2)*0.5 - 0.1260 - 0.3516*xp*2 + 0.2843*np.power(xp,2)*3 - 0.1015*np.power(xp,3)*4)
    return foil_der


delta_vec = [] # vector containing the computed deltas
# We iterate over file names
for x_fix in x:
    print('x_fix = {}'.format(x_fix))
    # File path to your CSV file
    csv_file = './{}/{}'.format(csv_folder.format(Re), csv_name.format(x_fix, Re))
    if slip:
        csv_file_slip = './{}/{}'.format(slip_csv_folder.format(Re), slip_csv_name.format(x_fix, Re))
        # Extraction of the velocity at the foil boundary from the slip simulation
        with open(csv_file_slip, 'r', newline='') as file_slip:
            reader = csv.DictReader(file_slip)
            u_slip_x = [] # x-component of the velocity in the slip case
            u_slip_y = [] # y-component of the velocity in the slip case
            for row in reader:
                u_slip_x.append(float(row['U:0']))
                u_slip_y.append(float(row['U:2']))
            
            # Convert to numpy array
            u_slip = np.array([u_slip_x, u_slip_y]) # vector of velocities along the normal in the slip case

    # Lists to store the extracted columns from the non-slip case
    U = [] # x-component to the velocity
    V = [] # y-component to the velocity
    X = [] # x-component of the points on the extracted line
    Y = [] # y-component of the points on the extracted line
    arcs = [] # arc length of the line 

    # Extraction of the velocity and points from the no-slip case
    with open(csv_file, 'r', newline='') as file:
        reader = csv.DictReader(file)

        for row in reader:
            U.append(float(row['U:0']))
            V.append(float(row['U:2']))
            X.append(float(row['Points:0']))
            Y.append(float(row['Points:2']))
            arcs.append(float(row['arc_length']))

    # Convert lists to NumPy arrays
    u = np.array([U, V]) # u[:,j] is the vector velocity at the j-th point
    points = np.array([X, Y]) # points[:,j] is the j-th point on the line
    arc_lengths = np.array(arcs) # arc_lengths[j] is the distance between points[:,j] and points[:,0]

    # Compute the tangent to the foil
    T = np.array([1, foil_derivative(t,chord,x_fix)])
    N = np.array([-foil_derivative(t,chord,x_fix), 1])
    T = T/np.linalg.norm(T) # unit vector tangent to the foil 
    N = N/np.linalg.norm(N) # unit normal vector to the foil
    print('  N = {}'.format(N))
    print('  T = {}'.format(T))

    # compute the parallel component of the velocity 
    u_par = np.dot(T,u)
    
    # compute U_e at the foil boundary if needed
    if slip:
        U_e = np.dot(u_slip[:,0], T) # parallel component of u_slip at the foil boundary
        print('  U_e = {}'.format(U_e))
        boundary_layer = Bl_thickness * U_e # boundary layer thickness in the slip case
    # otherwise compute u_max
    else:
        u_max = np.max(u_par) # u_max is defined as the maximum velocity tangent to the foil
        print('  u_max = {}'.format(u_max))
        boundary_layer = Bl_thickness * u_max # boundary layer thickness in the noslip case


    
    delta = 0.0 # initialize delta
    for j in range(len(u_par)):
        # When we reach the end of the boundary layer we compute the euclidian distance between 
        # A = points[:,0] and C = points[:,j]
        if u_par[j] >= boundary_layer: 
            # delta = np.dot(points[:,j] - points[:,0], N) # euclidian distance
            delta = arc_lengths[j]
            print('  A = {}'.format(points[:,0]))
            print('  C = {}'.format(points[:,j]))
            print('  delta = {}'.format(delta))
            break
    
    # We append the value of delta and x_fix to the arrays defined outside the loop over x_fix
    delta_vec.append(delta)

# Convert x and delta_vec into np.arrays
x = np.array(x)
delta_vec = np.array(delta_vec)

# plot the deltas
plt.scatter(x, delta_vec, color='red', zorder=5, label='Sample points')

# best fitting square root:
# Transform x -> sqrt(x)
sq = np.sqrt(x)
# Build design matrix: [u, 1]
A = np.vstack([sq, np.ones_like(sq)]).T
# Solve least squares for coefficients [a, b]
a, b = np.linalg.lstsq(A, delta_vec, rcond=None)[0]
# plot the fitting curve
plt.plot(np.linspace(x[0], x[-1], 100), a*np.sqrt(np.linspace(x[0], x[-1], 100)) +b, label='Fit: {:.5f}√x + {:.5f}'.format(a, b))

# theoretical curve
def f(x):
    return 4.9 * np.sqrt((chord/Re) * x)  # theoretical curve 
plt.plot(np.linspace(x[0], x[-1], 100), f(np.linspace(x[0], x[-1], 100)), label='Th: {:.5f}√x'.format(4.9 * np.sqrt(chord/Re)))


# Labels and legend
plt.xlabel('x')
plt.ylabel(r'$\delta(x)$')
if slip:
    plt.title('Boundary layer with slip simulation')
else:
    plt.title('Boundary layer without slip simulation')
plt.grid(True)
plt.legend()
plt.show()