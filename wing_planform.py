import matplotlib.pyplot as plt
from scipy.integrate import quad    
import math

# Adjust font size
plt.rcParams.update({'font.size': 22})

# --------------------------------------------- Parameters ---------------------------------------------

# Parameters
AR = 9                                          # [-]
S = 66.6                                        # [m^2]
c4_sweep = 26.3                                 # [deg]

# Convert parameters \ calculate if necessary
span = math.sqrt(AR*S)                          # [m]
taper_ratio = 0.2*(2-c4_sweep*(math.pi/180))    # [-] 
c_r = (2*S)/((1+taper_ratio)*span)              # [m]
c_t = taper_ratio*c_r                           # [m]
c4_sweep /= (180/math.pi)                       # [rad]

# -------------------------------------------- Calculations --------------------------------------------

# Calculation MAC: chord and location
# Integrate from root to tip 
ymin = 0
ymax = span/2

# Define function to integrate
def function(y): 
    return (c_r - ((c_r-c_t)/(span/2))*y)**2

# Integrate 
result, error = quad(function, ymin, ymax)

# Apply formulas to find MAC and spanwise position
mac = (2/S) * result
spanwise_pos = (mac - c_r) / -((c_r-c_t)/(span/2))

# Root chord position
q_c_r = 0
le_c_r = q_c_r + 0.25 * c_r
te_c_r = q_c_r - 0.75 * c_r

# Tip chord position
q_c_t = q_c_r - (span/2) * math.tan(c4_sweep)
le_c_t = q_c_t + 0.25 * c_t
te_c_t = q_c_t - 0.75 * c_t

# MAC position
q_c_mac = q_c_r - (spanwise_pos) * math.tan(c4_sweep)
le_mac = q_c_mac + 0.25 * mac
te_mac = q_c_mac - 0.75 * mac

# Calculating XLEMAC
XLEMAC = mac * math.tan(c4_sweep) + (1/4)*c_r - (1/4)*mac

# Calculate leading and trailing edge sweep angle
slope_le = (le_c_r - le_c_t) / (span/2)
le_sweep = math.atan(slope_le)

slope_te = (te_c_r - te_c_t) / (span/2)
te_sweep = math.atan(slope_te)

# ------------------------------------------ Printing results ------------------------------------------

print("-------------------- Summary ----------------------")
print("Span:                        ", round(span,2), "[m]")
print("Taper ratio:                 ", round(taper_ratio,3), "[-]")
print("Root chord:                  ", round(c_r,2), "[m]")
print("Tip chord:                   ", round(c_t,2), "[m]")
print("MAC:                         ", round(mac,2), "[m]")
print("Spanwise position MAC:       ", round(spanwise_pos,2), "[m]")
print("XLEMAC:                      ", round(XLEMAC,2), "[m]")
print("Leading edge sweep angle:    ", round(le_sweep*(180/math.pi),2), "[deg]")
print("Trailing edge sweep angle:   ", round(te_sweep*(180/math.pi),2), "[deg]")
print("---------------------------------------------------")

# ---------------------------------------------- Plotting ----------------------------------------------

# Make plot
plt.figure()

# Quarter chord line left side
plt.plot((0,-span/2),(q_c_r,q_c_t), label="Quarter Chord line", linewidth=3)

# Root chord
plt.plot((0,0),(le_c_r,te_c_r), label="Root chord", linewidth=3)

# Tip chords
plt.plot((-span/2,-span/2),(le_c_t,te_c_t), 'r', label="Tip chord", linewidth=3)
plt.plot((span/2,span/2),(le_c_t,te_c_t), 'r', linewidth=3)

# MAC on left side
mac_plot = plt.plot((-spanwise_pos,-spanwise_pos),(te_mac,le_mac), label="MAC", linewidth=3)

# Connecting
plt.plot((-span/2,0),(le_c_t,le_c_r), 'b', linewidth=3)
plt.plot((0,span/2),(le_c_r,le_c_t), 'b', linewidth=3)
plt.plot((-span/2,0),(te_c_t,te_c_r), 'b', linewidth=3)
plt.plot((0,span/2),(te_c_r,te_c_t), 'b', linewidth=3)

# Finalizing
plt.gca().set_aspect('equal', adjustable='box')
plt.title("Wing planform")
plt.xlabel("Spanwise position")
plt.ylabel("Longitudinal position")
plt.legend()
plt.grid()
plt.show()