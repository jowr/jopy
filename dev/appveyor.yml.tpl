
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
  - cmd: C:\Miniconda.exe /S /D=C:\Miniconda
  #- cmd: C:\Miniconda\Scripts\activate.bat root
  - cmd: SET PATH=C:\Miniconda3\Scripts;C:\Miniconda\Scripts;%PATH%
  # We need to do this first as other commands may not work with older versions of conda.
  - cmd: conda update -yq conda
  - cmd: conda install -yq {{ bas_pkgs }}
  - cmd: conda create -yq -n condaenv python=%PYTHON%

install:
  - cmd: activate condaenv
  - cmd: conda install -yq unxutils {{ cus_pkgs }}
  #- cmd: conda install -yq -c coolprop coolprop
  - cmd: pip install {{ pip_pkgs }}
  - cmd: python setup.py install

test_script:
  - cmd: nosetests --with-coverage --cover-package=jopy --logging-level=INFO
