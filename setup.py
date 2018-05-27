from setuptools import setup, find_packages

setup(
	name = 'spada',
	packages = find_packages(),
	package_dir = {'spada': 'spada'},
	package_data = {'spada': ['data/*.pkl']},
	version = '1.16',
	description = 'Find splicing-led, functional changes of the proteome. ',
	author = 'Héctor Climente-González',
	author_email = 'hector.climente@curie.fr',
	license = 'MIT',
	url = 'https://github.com/hclimente/spada',
	download_url = 'https://github.com/hclimente/spada/archive/v1.9.1.tar.gz',
	keywords = ['alternative', 'splicing', 'analysis', 'transcriptomics', 'networks'],
	classifiers = [
		'Development Status :: 4 - Beta',
		'Topic :: Scientific/Engineering :: Bio-Informatics',
		'Intended Audience :: Science/Research',
		'Intended Audience :: Healthcare Industry',
		'License :: OSI Approved :: MIT License',
		'Operating System :: POSIX :: Linux',
		'Operating System :: MacOS',
		'Operating System :: Microsoft :: Windows',
		'Programming Language :: Python :: 3 :: Only'],
	install_requires = [
		'networkx >= 1.11',
		'numpy >= 1.13.1'],
	scripts=['bin/spada']
)
