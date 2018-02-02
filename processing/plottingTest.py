
import DCIPtools as DCIP
import matplotlib.pyplot as plt

################################################################
# define the file required for import
fileName = "/Users/juan/Documents/testData/FieldSchool_2017new.DAT"

patch = DCIP.loadDias(fileName)   # Create the patch from data file

rdg = 0                           # Source to plot

f, plts = plt.subplots(2, 3)
plts[0, 0].plot(patch.window_center,
                patch.readings[rdg].Vdp[0].Vs /
                (patch.readings[rdg].Vdp[0].Vp / 1000.), 'b-o')
plts[0, 0].set_title("dipole 1")
plts[0, 1].plot(patch.window_center,
                patch.readings[rdg].Vdp[1].Vs /
                (patch.readings[rdg].Vdp[1].Vp / 1000.), 'b-o')
plts[0, 1].set_title("dipole 2")
plts[0, 2].plot(patch.window_center,
                patch.readings[rdg].Vdp[2].Vs /
                (patch.readings[rdg].Vdp[2].Vp / 1000.), 'b-o')
plts[0, 2].set_title("dipole 3")
plts[1, 0].plot(patch.window_center,
                patch.readings[rdg].Vdp[3].Vs /
                (patch.readings[rdg].Vdp[3].Vp / 1000.), 'b-o')
plts[1, 0].set_title("dipole 4")
plts[1, 1].plot(patch.window_center,
                patch.readings[rdg].Vdp[4].Vs /
                (patch.readings[rdg].Vdp[4].Vp / 1000.), 'b-o')
plts[1, 1].set_title("dipole 5")
plts[1, 2].plot(patch.window_center,
                patch.readings[rdg].Vdp[5].Vs /
                (patch.readings[rdg].Vdp[5].Vp / 1000.), 'b-o')
plts[1, 2].set_title("dipole 6")

plt.show()
