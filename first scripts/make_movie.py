import tqdm
import subprocess
import multiprocessing
import gc
import simulation
import matplotlib.pyplot as plt

plt.style.use('dark_background')

sim = simulation.MoriaSim('69002_p200.0_a600.0_r600.0_c8.15_movie', kicked=True)

sim.compute_cog(force=False, save_cache=True, use_process=True)
print(sim)

def save_fig(i):
    fig = sim.plot_gas(i, vmax=2e-1, vmin=1e-4, width=30, resolution=1000, sfh=True, cog=True)
    fig.savefig("{:03}.png".format(i), bbox_inches='tight', dpi=300)
    plt.close(fig)
    # del sim.snap_list[i]
    # gc.collect()

with multiprocessing.Pool(processes=7) as pool:

    print("Making snapshots:")
    # pool.map(save_fig, range(len(sim)))

    max_ = len(sim)
    with tqdm.tqdm(total=max_) as pbar:
        for i, _ in tqdm.tqdm(enumerate(pool.imap_unordered(save_fig, range(max_)))):
            pbar.update()

# for i in tqdm.tqdm(range(len(sim))): 
#     # print("saving ")
#     save_fig(i)

# for i in [0]:#tqdm.tqdm(range(len(sim))): 
#   # print("saving ")
#   fig = sim.plot_gas_and_stars(i, vmax=2e-1, vmin=5e-4, resolution=1000, sfh=True, cog=True)
#   fig.savefig("{:03}.png".format(i), bbox_inches='tight', dpi=300)
#   plt.close(fig)

scale = ( (4545//2)*2, (3257//2)*2 )
command = "mencoder mf://*.png -mf fps=10:type=png \
-ovc lavc -lavcopts vcodec=mpeg4:mbd=2:trell -vf scale={}:{} -oac copy -o output.mp4".format(*scale)

# Auto scale
# command = "mencoder mf://*.png -mf fps=10:type=png -ovc lavc -lavcopts vcodec=mpeg4:mbd=2:trell -v -oac copy -o output.mp4"

# ffmpeg
# command = 'ffmpeg -r 10 -f image2 -s 1920x1080 -i %03d.png -vcodec libx264 -crf 25 -pix_fmt yuv420p -vf "scale=trunc(iw/2)*2:trunc(ih/2)*2" test.mp4'

print(command)
subprocess.call(command, shell=True)
