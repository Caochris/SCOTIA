from setuptools import setup
import os

os.environ['BEZIER_NO_EXTENSION'] = 'True'

setup(
    name='dbscan-ot',
    version='1.0',
    packages=['dbscan_ot'],
    package_dir={'dbscan_ot': 'dbscan_ot'},
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
