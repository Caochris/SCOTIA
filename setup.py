from setuptools import setup
import os

os.environ['BEZIER_NO_EXTENSION'] = 'True'

setup(
    name='dbscan-ot',
    version='1.0',
    packages=['dbscan_ot'],
    package_dir={'dbscan_ot': 'dbscan_ot'},
    install_requires=[
        'bezier>=2020.5.19',
        'forest-fire-clustering>=0.0.25',
        'matplotlib>=3.3.4',
        'numpy',
        'pandas',
        'POT>=0.8.0',
        'seaborn',
        'scipy',
        'scanpy',
        'statsmodels',
        'scikit-learn'
    ]
)
