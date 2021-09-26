import io
import os
from setuptools import find_packages, setup


NAME = 'spinstats'
DESCRIPTION = 'Calculations of spin-statistical factors for triplet-triplet annihilation upconversion.'
URL = 'https://github.com/davidbossanyi/spin-statistics'
AUTHOR = 'David Bossanyi'
EMAIL = 'davebossanyi@gmail.com'
REQUIRES_PYTHON = '>=3.8'
VERSION = '0.0.1'

REQUIRED = [
    'numpy >=1.20', 'matplotlib >=3.3', 'scipy >=1.7', 'pandas >=1.3',
]

here = os.path.abspath(os.path.dirname(__file__))

try:
    with io.open(os.path.join(here, 'README.md'), encoding='utf-8') as f:
        long_description = '\n' + f.read()
except FileNotFoundError:
    long_description = DESCRIPTION

about = {}
if not VERSION:
    with open(os.path.join(here, NAME, '__version__.py')) as f:
        exec(f.read(), about)
else:
    about['__version__'] = VERSION

setup(
    name=NAME,
    version=about['__version__'],
    description=DESCRIPTION,
    long_description=long_description,
    long_description_content_type='text/markdown',
    author=AUTHOR,
    author_email=EMAIL,
    python_requires=REQUIRES_PYTHON,
    url=URL,
    packages=find_packages(exclude=('doc', 'examples', 'notebooks')),
    install_requires=REQUIRED,
    include_package_data=True,
    license='MIT',
    classifiers=[
        'License :: OSI Approved :: MIT License',
        'Intended Audience :: Science/Research',
        'Programming Language :: Python :: 3 :: Only',
        'Programming Language :: Python :: 3.6',
        'Operating System :: OS Independent',
        'Development Status :: 4 - Beta',
    ],
)
