from setuptools import setup, find_packages

setup(
  name='spada',
  packages=find_packages('.'),
  version='1.0',
  description='Find splicing-led, functional changes of the proteome. ',
  author='Héctor Climente-González',
  author_email='eduardo.eyras@upf.edu',
  license='MIT',
  url='https://github.com/hclimente/spada',
  download_url='https://github.com/hclimente/spada/archive/v1.0.tar.gz',
  keywords=['alternative', 'splicing', 'analysis', 'transcriptomics', 'networks'],
  classifiers=[
        'Development Status :: 5 - Production/Stable',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Operating System :: POSIX :: Linux',
        'Programming Language :: Python :: 3.5'],
  install_requires=['networkx>=1.11',
                    'pandas>=0.20.3',
					'numpy>=1.13.1']
)
