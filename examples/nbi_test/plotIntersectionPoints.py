import re
import numpy as np
import matplotlib.pyplot as plt

# --- Load and parse file ---
path = "beam_map_grouped.csv"

records = []
current_ir = None
ir_pat = re.compile(r"^\s*ir\s*=\s*(\d+)")
pt_pat = re.compile(
    r"s=([+-]?\d+(?:\.\d+)?(?:e[+-]?\d+)?)\s*,\s*r=([+-]?\d+(?:\.\d+)?(?:e[+-]?\d+)?)\s*,\s*theta=([+-]?\d+(?:\.\d+)?(?:e[+-]?\d+)?)",
    re.IGNORECASE,
)

with open(path, "r") as f:
    for line in f:
        m_ir = ir_pat.match(line)
        if m_ir:
            current_ir = int(m_ir.group(1))
            continue
        m_pt = pt_pat.search(line)
        if m_pt and current_ir is not None:
            s = float(m_pt.group(1))
            r = float(m_pt.group(2))
            theta = float(m_pt.group(3))
            records.append((current_ir, s, r, theta))

# Convert to numpy array: columns [ir, s, r, theta]
data = np.array(records)


# Find all unique flux surface indices:
def ir_choice(ir_choice):

    mask = data[:,0] == ir_choice

    s = data[mask,1].astype(float)
    r = data[mask,2].astype(float)
    theta = data[mask,3].astype(float)


    n  = np.array([0.0, 1.0, 0.0]) 
    P0 = np.array([0.685, -1.028, 0.0]) 
    

    e1 = np.array([0.0,  0.0,  1.0])
    e2 = np.array([1.0,   0.0,  0.0])

    # compute beam-local x/y offsets
    x_beam = r * np.cos(theta)
    y_beam = r * np.sin(theta)

    # vectorized conversion to Cartesian: P0 + s*n + e1*x_beam + e2*y_beam
    # shape: (N,3), i th
    ##i thino something is wrong here because z is not 0 when r is 0. should be able to do a simpler case
    points = (P0[np.newaxis, :] +
            s[:, np.newaxis] * n[np.newaxis, :] +
            x_beam[:, np.newaxis] * e1[np.newaxis, :] +
            y_beam[:, np.newaxis] * e2[np.newaxis, :])

    x = points[:, 0]
    y = points[:, 1]
    z = points[:, 2]
    
    tol = 0.001
    z_new = np.zeros_like(z)
    y_new = np.zeros_like(y)
    x_new = np.zeros_like(x)
    for i in range(len(z)):
        #project onto z=0 plane
        if z[i] < tol and z[i] > -tol:
            z_new[i] = z[i]
            x_new[i] = x[i]
            y_new[i] = y[i]

    radius_vector = np.linspace(0, 0.23-0.023, 10)
    #print(radius_vector)
    R0 = 0.798
    NR =len(radius_vector)
        
    r_f = radius_vector
    r_mid = 0.5 * (r_f[:-1] + r_f[1:])
    theta_new = np.linspace(-np.pi, np.pi, 100)
    colors = plt.cm.rainbow(np.linspace(0, 1, NR + 1))

    for i, r in enumerate(r_f):
        R1 = R0 + r
        R2 = R0 - r
        X1 = R1 * np.cos(theta_new)
        Y1 = R1 * np.sin(theta_new)
        X2 = R2 * np.cos(theta_new)
        Y2 = R2 * np.sin(theta_new)
        plt.plot(X1, Y1, '-', color=colors[i], label=f'ir = {i}')
        plt.plot(X2, Y2, '-', color=colors[i])
    
    plt.scatter(x_new, y_new, color=colors[ir_choice] , s= 5, label=f'Beam points ir={ir_choice}')
    plt.grid(True)
    plt.tight_layout()
    

for i in [1,2,3,4,5,6,7,8,9]:  #len(np.unique(data[:,0])):
    ir_choice(i)
    #plt.show()
#plt.legend()    
plt.show()


