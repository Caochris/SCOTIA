from setuptools import setup
import os

os.environ['BEZIER_NO_EXTENSION'] = 'True'

setup(
    name='scotia',
    version='1.0',
    packages=['scotia'],
    package_dir={'scotia': 'scotia'},
    license="MIT",
    description='scotia',
    url='https://github.com/Caochris/SCOTIA',
    python_requires='>=3.6',
    author='Jingyi Cao',
    author_email='jcao13@bwh.harvard.edu',
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
        'llvmlite',
        'xdoctest'
    ]
)
