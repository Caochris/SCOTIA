from setuptools import setup
import os

os.environ['BEZIER_NO_EXTENSION'] = 'True'

setup(
    name='scotia',
    version='1.0',
    packages=['scotia'],
    package_dir={'scotia': 'scotia'},
    install_requires=[
        'bezier',
        'forest-fire-clustering',
        'matplotlib',
        'numpy',
        'pandas',
        'POT',
        'seaborn',
        'scipy',
        'scanpy',
        'statsmodels',
        'scikit-learn',
        'llvmlite'
    ]
)
