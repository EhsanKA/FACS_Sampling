Python version 3.10

import pathlib
from setuptools import setup, find_packages

HERE = pathlib.Path(__file__).parent

VERSION = '0.1.0'
PACKAGE_NAME = 'FACS_Sampling'
AUTHOR = 'Ehsan Karimiara'
AUTHOR_EMAIL = 'ehsan.karimiara@mdc-berlin.de'
URL = 'https://github.com/EhsanKA/FACS_Sampling'

LICENSE = 'MIT License'
DESCRIPTION = 'New strategies on Sampling from FACS data with more focus on rare cell types.'

LONG_DESCRIPTION = (HERE / "README.md").read_text()
LONG_DESC_TYPE = "text/markdown"

INSTALL_REQUIRES = [
      'numpy',
      'pandas',
      'scanpy',
      'scipy',
]

setup(name=PACKAGE_NAME,
      version=VERSION,
      description=DESCRIPTION,
      long_description=LONG_DESCRIPTION,
      long_description_content_type=LONG_DESC_TYPE,
      author=AUTHOR,
      license=LICENSE,
      author_email=AUTHOR_EMAIL,
      url=URL,
      install_requires=INSTALL_REQUIRES,
      packages=find_packages()
      )