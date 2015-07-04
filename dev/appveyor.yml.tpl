
version: {{ version }}.{build}

# Tell appveyor to not use msbuild
build: false

environment:
  matrix:
    - PYTHON: 2.7
    - PYTHON: 3.3
    - PYTHON: 3.4

platform:
  - x86
  - x64

init:
  - "ECHO %PYTHON%"
  - ps: if($env:PYTHON -eq "2.7" -and $env:Platform -eq "x64"){Start-FileDownload 'http://repo.continuum.io/miniconda/Miniconda-latest-Windows-x86_64.exe' C:\Miniconda.exe; echo "Done"}
  - ps: if($env:PYTHON -eq "2.7" -and $env:Platform -eq "x86"){Start-FileDownload 'http://repo.continuum.io/miniconda/Miniconda-latest-Windows-x86.exe' C:\Miniconda.exe; echo "Done"}
  - ps: if($env:PYTHON -eq "3.3" -and $env:Platform -eq "x64"){Start-FileDownload 'http://repo.continuum.io/miniconda/Miniconda3-latest-Windows-x86_64.exe' C:\Miniconda.exe; echo "Done"}
  - ps: if($env:PYTHON -eq "3.3" -and $env:Platform -eq "x86"){Start-FileDownload 'http://repo.continuum.io/miniconda/Miniconda3-latest-Windows-x86.exe' C:\Miniconda.exe; echo "Done"}
  - ps: if($env:PYTHON -eq "3.4" -and $env:Platform -eq "x64"){Start-FileDownload 'http://repo.continuum.io/miniconda/Miniconda3-latest-Windows-x86_64.exe' C:\Miniconda.exe; echo "Done"}
  - ps: if($env:PYTHON -eq "3.4" -and $env:Platform -eq "x86"){Start-FileDownload 'http://repo.continuum.io/miniconda/Miniconda3-latest-Windows-x86.exe' C:\Miniconda.exe; echo "Done"}
  - C:\Miniconda.exe /S /D=C:\Miniconda
  - C:\Miniconda\Scripts\activate.bat root
  # We need to do this first as other commands may not work with older versions of conda.
  - "conda update -yq conda"
  - "conda install -yq {{ bas_pkgs }}"
  - "conda create -yq -n condaenv python=%PYTHON%"

install:
  - "activate condaenv"
  - "conda install -yq {{ cus_pkgs }}"
  #- "conda install -yq -c coolprop coolprop"
  - "pip install {{ pip_pkgs }}"
  - "python setup.py install"

test_script:
  - "nosetests --with-coverage --cover-package=jopy --logging-level=INFO"
