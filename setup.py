from setuptools import setup

setup(name='simulation',
      version='0.1',
      description='A module to analyze Gadget2 simulations',
   	  license="MIT",
   	  author='Michele Mastropietro',
   	  author_email='michele.mastropietro@ugent.be',
      install_requires=['pynbody', 'astropy'],
      py_modules=['simulation', 'glass_ic', 'sim_duration'],
      scripts=['scripts/get_maxid', 'scripts/read_header'],
      entry_points = {
        	'console_scripts': ['sim_duration=sim_duration:main',
                               'plot_sfh=plot_sfh:main',
                               'plot_trace=plot_trace:main']
    	}
      )