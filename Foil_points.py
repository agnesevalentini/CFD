# script to plot the normals to the foil given the points on the foil upper surface
import numpy as np
import matplotlib.pyplot as plt

# equation of the foil
def foil(t: float, chord:float, x):
    xp = x/chord
    return 5*t*chord*(0.2969*np.sqrt(xp) - 0.1260*xp - 0.3516*np.power(xp,2) + 0.2843*np.power(xp,3) - 0.1015*np.power(xp,4))

# derivative of the foil
def foil_derivative(t:float, chord:float, x):
    xp = x/chord
    foil_der = 5*t*(0.2969*np.power(xp,-1/2)*0.5 - 0.1260 - 0.3516*xp*2 + 0.2843*np.power(xp,2)*3 - 0.1015*np.power(xp,3)*4)
    return foil_der

r = 0.457452
h = 0.6
tau = np.linspace(0,1,200)
xaxis = np.linspace(0,0.8,200)

fig, ax = plt.subplots()
#inlet 
r_tau = 4*(h-r)*(tau-0.5)*(tau-0.5) + r
ax.plot(r_tau*np.cos(np.pi*(tau + 0.5)), r_tau*np.sin(np.pi*(tau + 1.5)), label='inlet')

#above and below
ax.plot(xaxis, np.full(200,h),'tab:blue')
ax.plot(xaxis, np.full(200,-h), 'tab:blue')

#foil
t = 0.12
chord = 0.48
x = np.linspace(0,chord,200)
ax.plot(x, foil(t, chord, x),'tab:green',label='foil')
ax.plot(x,-foil(t, chord, x),'tab:green')

# Tangente e Normale
X_Fix = [0.03, 0.06, 0.1, 0.2, 0.3, 0.4]
beta = 0.5
for x_fix in X_Fix:
    A = np.array([x_fix, foil(t, chord, x_fix)])
    T = np.array([1, foil_derivative(t,chord,x_fix)])
    N = np.array([-foil_derivative(t,chord,x_fix), 1])

    # final B, it's the one printed
    B = A + beta*N
    ax.plot([A[0], B[0]], [A[1], B[1]], 'tab:grey', linewidth=1)
    ax.plot([B[0]], [B[1]], 'tab:grey', marker='o', markersize=3)

    # intermediate B: one every 0.5
    betas = np.arange(0, beta-0.1, 0.1)
    for b in betas:
        B_i = A + b*N
        B_f = A + (b+0.1)*N
        ax.plot([B_i[0], B_f[0]], [B_i[1], B_f[1]], 'tab:red', marker='o', markersize=2, linewidth=1)

    print('x_fix = {}'.format(x_fix))
    print('  A = {}'.format(A))
    print('  B = {}'.format(B))
    #print('  t = {}'.format(T/np.linalg.norm(T)))
    #print('  n = {}'.format(N/np.linalg.norm(N)))



    ax.grid()
    ax.set_aspect('equal', adjustable='box')


plt.show()


