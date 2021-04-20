#!/usr/bin/env python3
import os
from setuptools import setup

directory = os.path.abspath(os.path.dirname(__file__))
with open(os.path.join(directory, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

setup(name='orbittools',
      version='0.1.0',
      description='Tools for orbit and mission analysis',
      author='Josh Stoffel',
      license='MIT',
      long_description=long_description,
      long_description_content_type='text/markdown',
      packages=['orbittools'],
      classifiers=[
          "Programming Language :: Python :: 3",
          "License :: OSI Approved :: MIT License"
      ],
      install_requires=['numpy', 'casadi'],
      python_requires='>=3.8')
