name='networkit'

version='5.0.3'

url='https://networkit.github.io/'

download_url='https://pypi.python.org/pypi/networkit'

license='MIT'

author='Christian L. Staudt, Henning Meyerhenke'

author_email = 'christian.staudt@kit.edu, meyerhenke@kit.edu'

description = 'NetworKit is a toolbox for high-performance network analysis'

long_description = """
NetworKit is a growing open-source toolkit for high-performance network analysis.
Its aim is to provide tools for the analysis of large networks in the size range
from thousands to billions of edges. For this purpose, it implements efficient
graph algorithms, many of them parallel to utilize multicore architectures. These
are meant to compute standard measures of network analysis, such as degree
sequences, clustering coefficients and centrality. In
this respect, NetworKit is comparable to packages such as NetworkX, albeit with a
focus on parallelism and scalability. NetworKit is also a testbed for algorithm
engineering and contains a few novel algorithms from recently published
research, especially in the area of community detection."""

keywords = ['graph algorithm', 'network analysis', 'social network']

platforms = 'any'

classifiers = [
'Development Status :: 5 - Production/Stable',
'Environment :: Console',
'Environment :: Other Environment',
'Framework :: IPython',
'Intended Audience :: Developers',
'Intended Audience :: End Users/Desktop',
'Intended Audience :: Science/Research',
'License :: OSI Approved :: MIT License',
'Natural Language :: English',
'Operating System :: OS Independent', 'Programming Language :: C++',
'Programming Language :: Python :: 3.3',
'Programming Language :: Python :: 3.4',
'Programming Language :: Python :: 3.5',
'Programming Language :: Python :: 3.6',
'Topic :: Software Development :: Libraries :: Python Modules',
'Topic :: Scientific/Engineering :: Bio-Informatics',
'Topic :: Scientific/Engineering :: Chemistry',
'Topic :: Scientific/Engineering :: Information Analysis',
'Topic :: Scientific/Engineering :: Mathematics',
]

install_requires = [
	'scipy',
	'matplotlib',
	'pandas',
	'numpy',
	'networkx',
	'tabulate',
	'ipython'
]
