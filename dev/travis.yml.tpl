# Config file for automatic testing at travis-ci.org

language: python

python:
  - "3.4"
  - "3.3"
  - "2.7"
  - "2.6"
  
#  - "pypy"

## command to install dependencies, e.g. pip install -r requirements.txt --use-mirrors
#install: 
#  - pip install -r requirements.txt
#
## command to run tests, e.g. python setup.py test
#script:
#  - py.test

notifications:
  email: false
 
# Update system and install miniconda, fill in the rest via pip
before_install:
  - sudo apt-get update -qq
  - sudo apt-get install -q -y gfortran libhdf5-serial-dev libnetcdf-dev liblapack-dev libatlas-dev
  - wget http://repo.continuum.io/miniconda/Miniconda-latest-Linux-x86_64.sh -O miniconda.sh
  - chmod +x miniconda.sh
  - ./miniconda.sh -b
  - export PATH=/home/travis/miniconda/bin:$PATH
  - conda update -yq conda
  - conda install -yq {{ bas_pkgs }}
  - conda create -yq -n condaenv python=$TRAVIS_PYTHON_VERSION
  # The next couple lines fix a crash with multiprocessing on Travis and are not specific to using Miniconda
  - sudo rm -rf /dev/shm
  - sudo ln -s /run/shm /dev/shm
  
# Install packages
install:
  - source activate condaenv
  - conda install -yq {{ cus_pkgs }}
  #- conda install -yq -c coolprop coolprop
  - pip install {{ pip_pkgs }}
  - python setup.py install
 
# Run test
script:
  - nosetests --with-coverage --cover-package=jopy --logging-level=INFO
 
# Calculate coverage
after_success:
  - coveralls
  - codecov

