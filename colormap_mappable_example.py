import numpy as np
import matplotlib
import matplotlib.pyplot as plt

n_plots = 10
how_many_colors = 12
cmap_name = 'jet'

norm = matplotlib.colors.Normalize(vmin=0, vmax=n_plots)
cmap = matplotlib.cm.get_cmap(cmap_name, how_many_colors)
mappable = plt.cm.ScalarMappable(norm=norm, cmap=cmap)

data_list = list()
for i in range(n_plots):
  data_list.append((np.random.rand(5), np.random.rand(5)))

fig, ax = plt.subplots()

for i in range(n_plots):
  ax.plot(*data_list[i], color=mappable.to_rgba(i))
fig.colorbar(mappable)
plt.show()