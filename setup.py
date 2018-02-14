from setuptools import setup

setup(name='simulation',
      version='0.1',
      description='A module to analyze Gadget2 simulations',
   	  license="MIT",
   	  author='Michele Mastropietro',
   	  author_email='michele.mastropietro@ugent.be',
      install_requires=['pynbody', 'astropy'],
      py_modules=['simulation'],
      )