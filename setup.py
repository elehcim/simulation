from setuptools import setup

setup(name='simulation',
      version='0.2',
      description='A package to analyze Gadget2 simulations',
      license="MIT",
      author='Michele Mastropietro',
      author_email='michele.mastropietro@ugent.be',
      install_requires=['pynbody', 'astropy', 'pandas', 'matplotlib', 'scipy', 'numpy',
                        'scikit-image', 'numpy-quaternion', 'numba', 'tqdm', 'photutils'],
      packages=['simulation', 'simulation.plot', 'simulation.parsers', 'simulation.interp'],
      package_data={'simulation': ['interp/tables/HI_data_all.npz', 'observationalData/*', 'data/*']},
      scripts=['simulation/scripts/get_maxid', 'simulation/scripts/read_header'],
      entry_points={
          'console_scripts': ['sim_duration=simulation.sim_duration:main',
                              'plot_sfh=simulation.plot.plot_sfh:main',
                              'plot_trace=simulation.plot.plot_trace:main',
                              'plot_dt=simulation.plot.plot_dt:main',
                              'plot_adhoc=simulation.plot.plot_adhoc:main',
                              'ssam=simulation.oo_lambda:main',
                              'omega_write=simulation.derotating_box:main',
                              'save_maps=simulation.save_maps:main',
                              'save_maps_all_bands=simulation.save_maps_all_bands:main',
                              'save_maps_hi=simulation.save_maps_hi:main',
                              'save_quat=simulation.save_quat:main',
                              'save_center=simulation.save_center:main',
                              'derotate_simulation=simulation.derotate_simulation:main',
                              'save_data=simulation.save_data:main',
                              'snap2mat=simulation.convert_snapshot_to_mat:main',
                              ]
                   }
      )
