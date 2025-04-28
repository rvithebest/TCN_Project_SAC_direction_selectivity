from neuron import h, gui
import matplotlib.pyplot as plt
import neuron
import numpy as np
BL_norm=-70
def entropy_calculator(p: np.ndarray) -> float:
    return -np.sum(np.where(p > 0, p * np.log(p), 0)) # to avoid log(0) issues
# h.load_file("new_simulation_excite.hoc")
h.load_file("SAC_inhibit_simulator_2.hoc")
# Direction: 1 : Right to left
current_direc = 1
neuron.hoc.execute(f"current_direc={current_direc}")
neuron.hoc.execute("update_direction()")
# Create a Vector to store the voltage
v_rec_soma = h.Vector()
# Right Quadrant distal dendrties- 92,93,99,116
v_rec_dend_92 = h.Vector()
v_rec_dend_93 = h.Vector()
v_rec_dend_99 = h.Vector()
v_rec_dend_116 = h.Vector()
# Left Quadrant distal dendrties- 261,262, 272, 273
v_rec_dend_261 = h.Vector()
v_rec_dend_262 = h.Vector()
v_rec_dend_272 = h.Vector()
v_rec_dend_273 = h.Vector()
t_rec = h.Vector()

# Record voltage from the soma
v_rec_soma.record(h.soma(0.5)._ref_v)  # Soma is assumed to be the main section
# Record voltage from distal dendrites
v_rec_dend_92.record(h.dend[92](0.5)._ref_v)
v_rec_dend_93.record(h.dend[93](0.5)._ref_v)
v_rec_dend_99.record(h.dend[99](0.5)._ref_v)
v_rec_dend_116.record(h.dend[116](0.5)._ref_v)
v_rec_dend_261.record(h.dend[261](0.5)._ref_v)
v_rec_dend_262.record(h.dend[262](0.5)._ref_v)
v_rec_dend_272.record(h.dend[272](0.5)._ref_v)
v_rec_dend_273.record(h.dend[273](0.5)._ref_v)
t_rec.record(h._ref_t)  # Record the time

# Set up the simulation parameters
h.tstop = 500  # Set simulation time (in ms)
h.dt = 0.025  # Set time step (in ms)
# Run the simulation
h.run()
# Right dendrites (subtract BL_norm from the max voltage for each)
dend_92_v_max = max(v_rec_dend_92) - BL_norm
dend_93_v_max = max(v_rec_dend_93) - BL_norm
dend_99_v_max = max(v_rec_dend_99) - BL_norm
dend_116_v_max = max(v_rec_dend_116) - BL_norm
# Left dendrites
dend_261_v_max = max(v_rec_dend_261) - BL_norm
dend_262_v_max = max(v_rec_dend_262) - BL_norm
dend_272_v_max = max(v_rec_dend_272) - BL_norm
dend_273_v_max = max(v_rec_dend_273) - BL_norm
# response- 0 or 1 (spike or no spike)
# 0 if v_rec is less than -10 or 1 if v_rec is greater than -10- depeding on the quadrant
# For Mutual information estimation
response_92_1 = 0 if (dend_92_v_max > 60) else 1
response_93_1 = 0 if (dend_93_v_max > 60) else 1
response_99_1 = 0 if (dend_99_v_max > 60) else 1
response_116_1 = 0 if (dend_116_v_max > 60) else 1
response_261_1 = 0 if (dend_261_v_max < 60) else 1
response_262_1 = 0 if (dend_262_v_max < 60) else 1
response_272_1 = 0 if (dend_272_v_max < 60) else 1
response_273_1 = 0 if (dend_273_v_max < 60) else 1
# Make two figures
# First figure- Simulatneous volatge trae of soma, dend 273 and dend 99
# with differet colors
plt.figure(figsize=(10, 5))
plt.plot(t_rec, v_rec_soma, label="Soma Voltage", color='black')
plt.plot(t_rec, v_rec_dend_261, label="Dend-Left Quadrant", color='red')
plt.plot(t_rec, v_rec_dend_99, label="Dend- Right Quadrant", color='blue')
plt.xlabel("Time (ms)", fontsize=15)
plt.ylabel("Membrane Potential (mV)", fontsize=15)
plt.title("Bar Moving from right to left- with inhibitory inputs", fontsize=16)
# plt.title("Bar Moving from right to left- Differential excitation (temporally)", fontsize=16)
# Make the font size of the title and labels bigger
plt.legend()
plt.grid()
plt.tight_layout()
plt.show()
# Save the figure-
# plt.savefig("./New_figures/fig_1_varicos.png", dpi=300, bbox_inches='tight')
# Second figure-8 subplots- corresponding to the 8 distal dendrites
plt.figure(figsize=(10, 5))
plt.subplot(4, 2, 1)
plt.plot(t_rec, v_rec_dend_92, label="Dend-Right Quadrant 92", color='blue')
plt.legend()
plt.ylim(-90, 40)
plt.subplot(4, 2, 2)
plt.plot(t_rec, v_rec_dend_261, label="Dend-Left Quadrant 261", color='red')
plt.legend()
plt.ylim(-90, 40)
plt.subplot(4, 2, 3)
plt.plot(t_rec, v_rec_dend_93, label="Dend-Right Quadrant 93", color='blue')
plt.legend()
plt.ylim(-90, 40)
plt.ylabel("Membrane Potential (mV)", fontsize=15)
plt.subplot(4, 2, 4)
plt.plot(t_rec, v_rec_dend_262, label="Dend-Left Quadrant 262", color='red')
plt.legend()
plt.ylim(-90, 40)
plt.subplot(4, 2, 5)
plt.plot(t_rec, v_rec_dend_99, label="Dend-Right Quadrant 99", color='blue')
plt.legend()
plt.ylim(-90, 40)
plt.subplot(4, 2, 6)
plt.plot(t_rec, v_rec_dend_272, label="Dend-Left Quadrant 272", color='red')
plt.legend()
plt.ylim(-90, 40)
plt.subplot(4, 2, 7)
plt.plot(t_rec, v_rec_dend_116, label="Dend-Right Quadrant 116", color='blue')
plt.xlabel("Time (ms)", fontsize=15)
plt.legend()
plt.ylim(-90, 40)
plt.subplot(4, 2, 8)
plt.plot(t_rec, v_rec_dend_273, label="Dend-Left Quadrant 273", color='red')
plt.xlabel("Time (ms)", fontsize=15)
plt.legend()
# All the y-axis range- shoubld be same- -70 to +20
plt.ylim(-90, 40)
plt.suptitle("Bar Moving from right to left with inhibitory inputs", fontsize=16)
# plt.suptitle("Bar Moving from right to left -Differential excitation (temporally)", fontsize=16)
# Kv3 conductance=0.01 S/cm^2
plt.show()
# Save the figure
# plt.savefig("./New_figures/fig_2_varicos.png", dpi=300)
# Direction: -1: Left to right
current_direc = -1
neuron.hoc.execute(f"current_direc={current_direc}")
neuron.hoc.execute("update_direction()")
# Create a Vector to store the voltage
v_rec_soma = h.Vector()
# Right Quadrant distal dendrties- 92,93,99,116
v_rec_dend_92 = h.Vector()
v_rec_dend_93 = h.Vector()
v_rec_dend_99 = h.Vector()
v_rec_dend_116 = h.Vector()
# Left Quadrant distal dendrties- 261,262, 272, 273
v_rec_dend_261 = h.Vector()
v_rec_dend_262 = h.Vector()
v_rec_dend_272 = h.Vector()
v_rec_dend_273 = h.Vector()
t_rec = h.Vector()

# Record voltage from the soma
v_rec_soma.record(h.soma(0.5)._ref_v)  # Soma is assumed to be the main section
# Record voltage from distal dendrites
v_rec_dend_92.record(h.dend[92](0.5)._ref_v)
v_rec_dend_93.record(h.dend[93](0.5)._ref_v)
v_rec_dend_99.record(h.dend[99](0.5)._ref_v)
v_rec_dend_116.record(h.dend[116](0.5)._ref_v)
v_rec_dend_261.record(h.dend[261](0.5)._ref_v)
v_rec_dend_262.record(h.dend[262](0.5)._ref_v)
v_rec_dend_272.record(h.dend[272](0.5)._ref_v)
v_rec_dend_273.record(h.dend[273](0.5)._ref_v)
t_rec.record(h._ref_t)  # Record the time

# Set up the simulation parameters
h.tstop = 500  # Set simulation time (in ms)
h.dt = 0.025  # Set time step (in ms)
# Run the simulation
h.run()
# Right dendrites (subtract BL_norm from the max voltage for each)
dend_92_v_max_2 = max(v_rec_dend_92) - BL_norm
dend_93_v_max_2 = max(v_rec_dend_93) - BL_norm
dend_99_v_max_2 = max(v_rec_dend_99) - BL_norm
dend_116_v_max_2 = max(v_rec_dend_116) - BL_norm
# Left dendrites
dend_261_v_max_2 = max(v_rec_dend_261) - BL_norm
dend_262_v_max_2 = max(v_rec_dend_262) - BL_norm
dend_272_v_max_2 = max(v_rec_dend_272) - BL_norm
dend_273_v_max_2 = max(v_rec_dend_273) - BL_norm
# Response- 0 or 1 (spike or no spike)
# 0 if v_rec is less than -10 or 1 if v_rec is greater than -10
# For Mutual information estimation
response_92_2 = 0 if (dend_92_v_max_2 > 60) else 1
response_93_2 = 0 if (dend_93_v_max_2 > 60) else 1
response_99_2 = 0 if (dend_99_v_max_2 > 60) else 1
response_116_2 = 0 if (dend_116_v_max_2 > 60) else 1
response_261_2 = 0 if (dend_261_v_max_2 < 60) else 1
response_262_2 = 0 if (dend_262_v_max_2 < 60) else 1
response_272_2 = 0 if (dend_272_v_max_2 < 60) else 1
response_273_2 = 0 if (dend_273_v_max_2 < 60) else 1
# Calculate the response for each dendrite
# Calculate Direction Selectivity Index (DSI) for each dendrite
DSI_92 = abs(dend_92_v_max - dend_92_v_max_2) / (dend_92_v_max + dend_92_v_max_2)
DSI_93 = abs(dend_93_v_max - dend_93_v_max_2) / (dend_93_v_max + dend_93_v_max_2)
DSI_99 = abs(dend_99_v_max - dend_99_v_max_2) / (dend_99_v_max + dend_99_v_max_2)
DSI_116 = abs(dend_116_v_max - dend_116_v_max_2) / (dend_116_v_max + dend_116_v_max_2)
DSI_261 = abs(dend_261_v_max - dend_261_v_max_2) / (dend_261_v_max + dend_261_v_max_2)
DSI_262 = abs(dend_262_v_max - dend_262_v_max_2) / (dend_262_v_max + dend_262_v_max_2)
DSI_272 = abs(dend_272_v_max - dend_272_v_max_2) / (dend_272_v_max + dend_272_v_max_2)
DSI_273 = abs(dend_273_v_max - dend_273_v_max_2) / (dend_273_v_max + dend_273_v_max_2)
# Create a list of all DSI values for the current Ra
DSI_values = [DSI_92, DSI_93, DSI_99, DSI_116, DSI_261, DSI_262, DSI_272, DSI_273]

# Calculate the mean and std of DSI values for this Ra
mean_DSI = np.mean(DSI_values)
std_DSI = np.std(DSI_values)
print(f"Mean DSI: {mean_DSI}, Std DSI: {std_DSI}")
# Mutual information
prob_1_given_1 = (response_92_1 + response_93_1 + response_99_1 + response_116_1 + response_261_1 + response_262_1 + response_272_1 + response_273_1) / 8
prob_1_given_2 = ( response_92_2 + response_93_2 + response_99_2 + response_116_2 + response_261_2 + response_262_2 + response_272_2 + response_273_2) / 8
prob_0_given_1 = 1 - prob_1_given_1
prob_0_given_2 = 1 - prob_1_given_2
prob_1 = (prob_1_given_1 + prob_1_given_2) / 2
prob_0 = 1 - prob_1
H_R = entropy_calculator(np.array([prob_1, prob_0]))
H_R_given_1 = entropy_calculator(np.array([prob_1_given_1, prob_0_given_1]))
H_R_given_2 = entropy_calculator(np.array([prob_1_given_2, prob_0_given_2]))
H_R_given_S = (H_R_given_1 + H_R_given_2) / 2  # Assuming P(S=1) = P(S=2) = 0.5
MI = H_R - H_R_given_S
print(f"Mutual Information: {MI}")
# Make two figures
# First figure- SImulatneous volatge trae of soma, dend 273 and dend 99
# with differet colors
plt.figure(figsize=(10, 5))
plt.plot(t_rec, v_rec_soma, label="Soma Voltage", color='black')
plt.plot(t_rec, v_rec_dend_273, label="Dend-Left Quadrant", color='red')
plt.plot(t_rec, v_rec_dend_99, label="Dend- Right Quadrant", color='blue')
plt.xlabel("Time (ms)", fontsize=15)
plt.ylabel("Membrane Potential (mV)", fontsize=15)
plt.title("Voltage Traces (Bar Moving from left to right)- with inhibitory inputs", fontsize=16)
# plt.title("Bar Moving from left to right- Differential excitation (temporally)", fontsize=16)
# Make the font size of the title and labels bigger
plt.legend()
plt.grid()
plt.tight_layout()
plt.show()
# Save the figure-
#plt.savefig("./New_figures/fig_3_varicos.png", dpi=300, bbox_inches='tight')
# Second figure-8 subplots- corresponding to the 8 distal dendrites
plt.figure(figsize=(10, 5))
plt.subplot(4, 2, 1)
plt.plot(t_rec, v_rec_dend_92, label="Dend-Right Quadrant 92", color='blue')
plt.legend()
plt.ylim(-90, 40)
plt.subplot(4, 2, 2)
plt.plot(t_rec, v_rec_dend_261, label="Dend-Left Quadrant 261", color='red')
plt.legend()
plt.ylim(-90, 40)
plt.subplot(4, 2, 3)
plt.plot(t_rec, v_rec_dend_93, label="Dend-Right Quadrant 93", color='blue')
plt.legend()
plt.ylim(-90, 40)
plt.ylabel("Membrane Potential (mV)", fontsize=15)
plt.subplot(4, 2, 4)
plt.plot(t_rec, v_rec_dend_262, label="Dend-Left Quadrant 262", color='red')
plt.legend()
plt.ylim(-90, 40)
plt.subplot(4, 2, 5)
plt.plot(t_rec, v_rec_dend_99, label="Dend-Right Quadrant 99", color='blue')
plt.legend()
plt.ylim(-90, 40)
plt.subplot(4, 2, 6)
plt.plot(t_rec, v_rec_dend_272, label="Dend-Left Quadrant 272", color='red')
plt.legend()
plt.ylim(-90, 40)
plt.subplot(4, 2, 7)
plt.plot(t_rec, v_rec_dend_116, label="Dend-Right Quadrant 116", color='blue')
plt.xlabel("Time (ms)", fontsize=15)
plt.legend()
plt.ylim(-90, 40)
plt.subplot(4, 2, 8)
plt.plot(t_rec, v_rec_dend_273, label="Dend-Left Quadrant 273", color='red')
plt.xlabel("Time (ms)", fontsize=15)
plt.legend()
# All the y-axis range- shoubld be same- -70 to +20
plt.ylim(-90, 40)
plt.suptitle("Bar Moving from left to right- with inhibitory inputs", fontsize=16)
# plt.suptitle("Bar Moving from left to right- Differential excitation (temporally)", fontsize=16)
plt.show()
# Save the figure
# plt.savefig("./New_figures/fig_4_varicos.png", dpi=300)


# Function to plot morphology

def plot_morphology(ax):
    for sec in h.allsec():
        n3d = int(h.n3d(sec=sec))
        if n3d > 0:
            xs = [h.x3d(i, sec=sec) for i in range(n3d)]
            ys = [h.y3d(i, sec=sec) for i in range(n3d)]
            ax.plot(xs, ys, color='black', linewidth=1)

# Function to get (x,y,z) at a given location on a section
def get_xyz_at(sec, loc):
    n3d = int(h.n3d(sec=sec))
    arc_lengths = [h.arc3d(i, sec=sec) for i in range(n3d)]
    total_length = arc_lengths[-1]
    target = loc * total_length
    for i in range(n3d - 1):
        if arc_lengths[i] <= target <= arc_lengths[i+1]:
            frac = (target - arc_lengths[i]) / (arc_lengths[i+1] - arc_lengths[i])
            x = h.x3d(i, sec=sec) + frac * (h.x3d(i+1, sec=sec) - h.x3d(i, sec=sec))
            y = h.y3d(i, sec=sec) + frac * (h.y3d(i+1, sec=sec) - h.y3d(i, sec=sec))
            z = h.z3d(i, sec=sec) + frac * (h.z3d(i+1, sec=sec) - h.z3d(i, sec=sec))
            return x, y, z
    return h.x3d(0, sec=sec), h.y3d(0, sec=sec), h.z3d(0, sec=sec)  # fallback

# --------
# Now extract all synapse locations automatically
synapse_positions = []
for i in range(int(h.synlist.count())):
    syn = h.synlist.o(i)
    sec = syn.get_segment().sec
    loc = syn.get_segment().x
    synapse_positions.append((sec, loc))

# --------
# Plot
fig = plt.figure(figsize=(8, 8))
ax = fig.add_subplot(111)

# Plot morphology
plot_morphology(ax)

# Plot synapses
for sec, loc in synapse_positions:
    x, y, z = get_xyz_at(sec, loc)
    ax.plot(x, y, 'ro', markersize=6)  # red circles for synapses

ax.set_aspect('equal')
ax.set_title('SAC Morphology with Synapse Locations', fontsize=16)
ax.set_xlabel('X (µm)', fontsize=14)
ax.set_ylabel('Y (µm)', fontsize=14)
plt.grid(True)
# plt.show()
# Save the figure
# plt.savefig("./New_figures/SAC_morphology_synapses_160.png", dpi=300, bbox_inches='tight')