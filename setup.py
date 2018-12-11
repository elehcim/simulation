from setuptools import setup, find_packages

# print(find_packages())

setup(name='simulation',
      version='0.1',
      description='A module to analyze Gadget2 simulations',
      license="MIT",
      author='Michele Mastropietro',
      author_email='michele.mastropietro@ugent.be',
      install_requires=['pynbody', 'astropy', 'pandas', 'matplotlib', 'scipy', 'numpy'],
      py_modules=['glass_ic', 'sim_duration', 'luminosity'],
      # packages=find_packages(),
      scripts=['scripts/get_maxid', 'scripts/read_header'],
      entry_points={
          'console_scripts': ['sim_duration=sim_duration:main',
                              'plot_sfh=plot.plot_sfh:main',
                              'plot_trace=plot.plot_trace:main']
      }
      )
