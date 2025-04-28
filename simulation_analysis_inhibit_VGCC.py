import numpy as np
import matplotlib.pyplot as plt
from neuron import h
import neuron

# Load the simulation file
# h.load_file("new_simulation_inhibit.hoc")
h.load_file("SAC_inhibit_simulator_2.hoc")
# Set up the simulation parameters
h.tstop = 500  # Set simulation time (in ms)
h.dt = 0.025  # Set time step (in ms)
VGCC_con=[0.015,0.02,0.025,0.03,0.035,0.04,0.045,0.05] # S/cm^2
# Lists to store the mean and std values for each Ra across dendrites
mean_DSI_values = []
std_DSI_values = []
mutual_info_values = []
BL_norm= -70  # Resting membrane potential in mV
def entropy_calculator(p: np.ndarray) -> float:
    return -np.sum(np.where(p > 0, p * np.log2(p), 0)) # to avoid log(0) issues
# Loop through Ra values
for current_con in VGCC_con:
    print(f"Running simulation for current_con = {current_con} S/cm^2")
    # Create Vectors to store voltage data
    v_rec_soma = h.Vector()
    v_rec_dend_92 = h.Vector()
    v_rec_dend_93 = h.Vector()
    v_rec_dend_99 = h.Vector()
    v_rec_dend_116 = h.Vector()
    v_rec_dend_261 = h.Vector()
    v_rec_dend_262 = h.Vector()
    v_rec_dend_272 = h.Vector()
    v_rec_dend_273 = h.Vector()
    t_rec = h.Vector()

    # Record voltage from soma and dendrites
    v_rec_soma.record(h.soma(0.9)._ref_v)
    v_rec_dend_92.record(h.dend[92](0.9)._ref_v)
    v_rec_dend_93.record(h.dend[93](0.9)._ref_v)
    v_rec_dend_99.record(h.dend[99](0.9)._ref_v)
    v_rec_dend_116.record(h.dend[116](0.9)._ref_v)
    v_rec_dend_261.record(h.dend[261](0.9)._ref_v)
    v_rec_dend_262.record(h.dend[262](0.9)._ref_v)
    v_rec_dend_272.record(h.dend[272](0.9)._ref_v)
    v_rec_dend_273.record(h.dend[273](0.9)._ref_v)
    t_rec.record(h._ref_t)  # Record time

    # Update the Ra value in the hoc file
    neuron.hoc.execute(f"gcan={current_con}")
    neuron.hoc.execute("update_gcan_con()")
    # Update the direction in the hoc file
    current_direc = 1
    print(f"current_direc={current_direc}")
    neuron.hoc.execute(f"current_direc={current_direc}")
    neuron.hoc.execute("update_direction()")

    # Run the simulation for forward direction
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
    # 0 if v_rec is less than -10 or 1 if v_rec is greater than -10- depending on the quadrant
    # For Mutual information estimation
    response_92_1= 0 if (dend_92_v_max > 60) else 1
    response_93_1= 0 if (dend_93_v_max > 60) else 1
    response_99_1= 0 if (dend_99_v_max > 60) else 1
    response_116_1= 0 if (dend_116_v_max > 60) else 1
    response_261_1= 0 if (dend_261_v_max < 60) else 1
    response_262_1= 0 if (dend_262_v_max < 60) else 1
    response_272_1= 0 if (dend_272_v_max < 60) else 1
    response_273_1= 0 if (dend_273_v_max < 60) else 1
    # opposite direction
    # Create Vectors to store voltage data
    v_rec_soma = h.Vector()
    v_rec_dend_92 = h.Vector()
    v_rec_dend_93 = h.Vector()
    v_rec_dend_99 = h.Vector()
    v_rec_dend_116 = h.Vector()
    v_rec_dend_261 = h.Vector()
    v_rec_dend_262 = h.Vector()
    v_rec_dend_272 = h.Vector()
    v_rec_dend_273 = h.Vector()
    t_rec = h.Vector()

    # Record voltage from soma and dendrites
    v_rec_soma.record(h.soma(0.9)._ref_v)
    v_rec_dend_92.record(h.dend[92](0.9)._ref_v)
    v_rec_dend_93.record(h.dend[93](0.9)._ref_v)
    v_rec_dend_99.record(h.dend[99](0.9)._ref_v)
    v_rec_dend_116.record(h.dend[116](0.9)._ref_v)
    v_rec_dend_261.record(h.dend[261](0.9)._ref_v)
    v_rec_dend_262.record(h.dend[262](0.9)._ref_v)
    v_rec_dend_272.record(h.dend[272](0.9)._ref_v)
    v_rec_dend_273.record(h.dend[273](0.9)._ref_v)
    t_rec.record(h._ref_t)  # Record time


    # Change direction and run simulation for backward direction
    current_direc = -1
    print(f"current_direc={current_direc}")
    neuron.hoc.execute(f"current_direc={current_direc}")
    neuron.hoc.execute("update_direction()")
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
    response_92_2= 0 if (dend_92_v_max_2 > 60) else 1
    response_93_2= 0 if (dend_93_v_max_2 > 60) else 1
    response_99_2= 0 if (dend_99_v_max_2 > 60) else 1
    response_116_2= 0 if (dend_116_v_max_2 > 60) else 1
    response_261_2= 0 if (dend_261_v_max_2 < 60) else 1
    response_262_2= 0 if (dend_262_v_max_2 < 60) else 1
    response_272_2= 0 if (dend_272_v_max_2 < 60) else 1
    response_273_2= 0 if (dend_273_v_max_2 < 60) else 1
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

    # Append the mean and std to their respective lists
    mean_DSI_values.append(mean_DSI)
    std_DSI_values.append(std_DSI)
    # Mutual information
    prob_1_given_1= (response_92_1 + response_93_1 + response_99_1 + response_116_1+ response_261_1 + response_262_1 + response_272_1 + response_273_1)/8
    prob_1_given_2= (response_92_2 + response_93_2 + response_99_2 + response_116_2+ response_261_2 + response_262_2 + response_272_2 + response_273_2)/8
    prob_0_given_1= 1 - prob_1_given_1
    prob_0_given_2= 1 - prob_1_given_2
    prob_1= (prob_1_given_1 + prob_1_given_2)/2
    prob_0= 1 - prob_1
    H_R= entropy_calculator(np.array([prob_1, prob_0]))
    H_R_given_1= entropy_calculator(np.array([prob_1_given_1, prob_0_given_1]))
    H_R_given_2= entropy_calculator(np.array([prob_1_given_2, prob_0_given_2]))
    H_R_given_S= (H_R_given_1 + H_R_given_2)/2 # Assuming P(S=1) = P(S=2) = 0.5
    MI= H_R - H_R_given_S
    mutual_info_values.append(MI)


# Plot the mean and std DSI values as error bars
plt.figure(figsize=(10, 6))

# Using plt.errorbar to plot the mean DSI with std as error bars
plt.errorbar(VGCC_con, mean_DSI_values, yerr=std_DSI_values, fmt='o', label='Mean DSI ± Std', capsize=5)

plt.xlabel('N-type Calcium Channel Conductance (S/cm^2)')
plt.ylabel('DSI')
plt.title('Mean DSI ± Standard Deviation for Different N-type Calcium channel conductance (With Inhibitory inputs)')
plt.grid(True)
plt.legend()
# Label the coordinates of point with maximum DSI and maximum MI
max_dsi_index = np.argmax(mean_DSI_values)
max_dsi_value = mean_DSI_values[max_dsi_index]
max_dsi_con = VGCC_con[max_dsi_index]
plt.annotate(f'Max DSI: {max_dsi_value:.2f} at N-type cond={max_dsi_con}',xy=(max_dsi_con, max_dsi_value),xytext=(max_dsi_con -0.01, max_dsi_value -0.05),arrowprops=dict(facecolor='black', arrowstyle='->'),fontsize=10)
plt.axhline(y=max_dsi_value, color='r', linestyle='--', label='Max DSI')
plt.legend()
plt.show()
# Plot the mutual information values
plt.figure(figsize=(10, 6))
plt.plot(VGCC_con, mutual_info_values, marker='o', label='Mutual Information', color='orange')
plt.xlabel('Ra Values (Ohm-cm)')
plt.ylabel('Mutual Information')
plt.title('Mutual Information for Different N-type Calcium channel conductance (With Inhibitory inputs)')
plt.grid(True)
plt.legend()
# Label the coordinates of point with maximum MI
max_mi_index = np.argmax(mutual_info_values)
max_mi_value = mutual_info_values[max_mi_index]
max_mi_con = VGCC_con[max_mi_index]

plt.annotate(f'Max MI: {max_mi_value:.2f} at N-type cond={max_mi_con}',xy=(max_mi_con, max_mi_value),xytext=(max_mi_con + 0.01, max_mi_value -0.05),arrowprops=dict(facecolor='black', arrowstyle='->'),fontsize=10)
plt.axhline(y=max_mi_value, color='r', linestyle='--', label='Max MI')
plt.legend()
plt.show()
