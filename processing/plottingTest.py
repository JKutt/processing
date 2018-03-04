
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
print(num_half_T)
trim_length = num_half_T * samples_half_t
xt = xt[0:int(trim_length)]
xt = np.asarray(xt)
# Ax = np.reshape(xt, (int(samples_half_t), int(num_half_T)), order='F')
# end read =====================================================

# stack time ===================================================
window = DCIP.createHanningWindow(num_half_T - 2)  # creates filter window
tHK = DCIP.filterKernal(filtershape=window)       # creates filter kernal
amp = DCIP.getFrequnceyResponse(xt)
plt.plot(np.log(amp))
stack = tHK * xt                                  # stack data
# plt.plot(-stack)
# plt.show()
start_vp = 50                      # start of Vp calculation (%)
end_vp = 90                        # end of Vp calculation (%)
Vp = DCIP.getPrimaryVoltage(start_vp,
                            end_vp,
                            stack)  # calculate the Vp
print(Vp)
# # end stacking ==================================================

# # Vs Windowing Time =============================================
# create time windows for Vs
timeFrom = [2040, 2060, 2080, 2120, 2160, 2200,
            2240, 2320, 2400,
            2480, 2560, 2640,
            2720, 2800, 2960,
            3120, 3280, 3440,
            3600, 3760]
timeTo = [2060, 2080, 2120, 2160, 2200, 2240,
          2320, 2400, 2480, 2560, 2640, 2720,
          2800, 2960, 3120, 3280, 3440,
          3600, 3760, 3920]
dkernal = DCIP.decayKernal(window_starts=np.asarray(timeFrom),
                           window_ends=np.asarray(timeTo),
                           window_weight=401,
                           window_overlap=45)  # creates decay kernal
decay = dkernal * stack                   # calculates the decay
# plt.plot(dkernal.getWindowCenters(), -decay / (Vp), 'bo-')
# plt.show()
# # end windowing ===================================================
t_ = dkernal.getWindowWidths()
w_ = np.ones(t_.size)
decay = -decay / (Vp)             # normalise decay by primary voltage
c, tau, M, error, vs = DCIP.getColeCole(decay, 0.45, 0.0025, 1.0, t_, w_)
print(c, tau, M, error)
# plt.plot(timet, vs, '-ko')
plt.show()
