
import DCIPtools as DCIP
import matplotlib.pyplot as plt
import numpy as np

################################################################
# Read binary time-series ======================================
xyz = open("/Users/juan/Documents/testData/ts2new.xyz")
t = []
xt = []
samples_half_t = 600.0
for line in xyz:
    x, y = line.split()
    t.append(float(x))
    xt.append(float(y))
xyz.close()
xt = np.asarray(xt)
length_signal = xt.size
num_half_T = np.round(xt.size / samples_half_t) - 2
trim_length = num_half_T * samples_half_t
xt = xt[0:int(trim_length)]
xt = np.asarray(xt)
# Ax = np.reshape(xt, (int(samples_half_t), int(num_half_T)), order='F')
# end read =====================================================

# stack time ===================================================
window = DCIP.createHanningWindow(39)  # creates filter kernal
# Polarity = np.power(-1, 0)
tHK = DCIP.filterKernal(filtershape=window)
# stack = xt * tHK
stk = np.matmul(Ax, tHK)
print(tHK.kernal)
# stk = np.matmul(Ax, bHK)                # produce the stack
# Vp = DCIP.getPrimaryVoltage(50, 90, stack)
# print(Vp)
# # end stacking ==================================================

# # Vs Windowing Time =============================================
# # create time windows for Vs
# timeFrom = [2040, 2060, 2080, 2120, 2160, 2200,
#             2240, 2320, 2400,
#             2480, 2560, 2640,
#             2720, 2800, 2960,
#             3120, 3280, 3440,
#             3600, 3760]
# timeTo = [2060, 2080, 2120, 2160, 2200, 2240,
#           2320, 2400, 2480, 2560, 2640, 2720,
#           2800, 2960, 3120, 3280, 3440,
#           3600, 3760, 3920]
# decay = DCIP.getWeightedVs(stack,
#                            timeFrom,
#                            timeTo,
#                            401)        # calculates decay data
# timet = (np.asarray(timeTo) +
#          np.asarray(timeFrom)) / 2.    # creates window center
# t_ = (np.asarray(timeTo) -
#       np.asarray(timeFrom))            # creates window centers
# plt.plot(timet, decay, 'bo-')
# # plt.show()
# # end windowing ===================================================
# decay = decay / (Vp)             # normalise decay by primary voltage
# c, tau, M, error, vs = DCIP.calcColeCole(decay, t_)
# print(c, tau, M, error)
# # plt.plot(timet, vs, '-o')
# plt.show()

# define the file required for import
# fileName = "/Users/juan/Documents/testData/L6050-err1.DAT"

# patch = DCIP.loadDias(fileName)   # Create the patch from data file

# rdg = 0                           # Source to plot
# dpnum = []
# for i in range(len(patch.readings[rdg].Vdp)):
#     dpnum.append(patch.readings[rdg].Vdp[i].dipole)
#     plt.plot(patch.window_center,
#              patch.readings[rdg].Vdp[i].Vs /
#              (patch.readings[rdg].Vdp[i].Vp / 1000.), '-o')
# plt.legend(dpnum)
# plt.title("All Decay")
# plt.xlabel("time (ms)")
# plt.ylabel("Voltage (mV)")
# plt.show()
