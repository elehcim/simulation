import chyplot
import enums
input_file  = '/home/michele/sim/sim67001/snapshot_0036'

print("\n-- Initializing hyplot")

reader = chyplot.CDataGadget()
writer = chyplot.CWriteGadget()

reader.setFilename(input_file)

data = reader.readFile()

data.rcom(True, enums.T_star, 0, 0, 0, True)
data.vcom(True, enums.T_star)
galaxy = chyplot.CGalaxy(data)

r_eff = galaxy.getHalfLightRadius()

print r_eff