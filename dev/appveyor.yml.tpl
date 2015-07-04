# Tell appveyor to not use msbuild
build: false

environment:
  matrix:
    - PYTHON: 2.7
      PYCOSAT: 0.6.0
    - PYTHON: 2.7
      PYCOSAT:
    - PYTHON: 3.4
      PYCOSAT: 0.6.0
    - PYTHON: 3.4
      PYCOSAT:

platform:
  - x86
  - x64

init:
  - "ECHO %PYTHON%"
  - "ECHO %PYCOSAT%"
  - ps: if($env:PYTHON -eq "2.7" -and $env:Platform -eq "x64" -and $env:PYCOSAT -ne "0.6.0"){Start-FileDownload 'http://repo.continuum.io/miniconda/Miniconda-3.5.5-Windows-x86_64.exe' C:\Miniconda.exe; echo "Done"}
  - ps: if($env:PYTHON -eq "2.7" -and $env:Platform -eq "x64" -and $env:PYCOSAT -eq "0.6.0"){Start-FileDownload 'http://repo.continuum.io/miniconda/Miniconda-3.4.2-Windows-x86_64.exe' C:\Miniconda.exe; echo "Done"}
  - ps: if($env:PYTHON -eq "2.7" -and $env:Platform -eq "x86" -and $env:PYCOSAT -ne "0.6.0"){Start-FileDownload 'http://repo.continuum.io/miniconda/Miniconda-3.5.5-Windows-x86.exe' C:\Miniconda.exe; echo "Done"}
  - ps: if($env:PYTHON -eq "2.7" -and $env:Platform -eq "x86" -and $env:PYCOSAT -eq "0.6.0"){Start-FileDownload 'http://repo.continuum.io/miniconda/Miniconda-3.4.2-Windows-x86.exe' C:\Miniconda.exe; echo "Done"}
  - ps: if($env:PYTHON -eq "3.4" -and $env:Platform -eq "x64" -and $env:PYCOSAT -ne "0.6.0"){Start-FileDownload 'http://repo.continuum.io/miniconda/Miniconda3-3.5.5-Windows-x86_64.exe' C:\Miniconda.exe; echo "Done"}
  - ps: if($env:PYTHON -eq "3.4" -and $env:Platform -eq "x64" -and $env:PYCOSAT -eq "0.6.0"){Start-FileDownload 'http://repo.continuum.io/miniconda/Miniconda3-3.4.2-Windows-x86_64.exe' C:\Miniconda.exe; echo "Done"}
  - ps: if($env:PYTHON -eq "3.4" -and $env:Platform -eq "x86" -and $env:PYCOSAT -ne "0.6.0"){Start-FileDownload 'http://repo.continuum.io/miniconda/Miniconda3-3.5.5-Windows-x86.exe' C:\Miniconda.exe; echo "Done"}
  - ps: if($env:PYTHON -eq "3.4" -and $env:Platform -eq "x86" -and $env:PYCOSAT -eq "0.6.0"){Start-FileDownload 'http://repo.continuum.io/miniconda/Miniconda3-3.4.2-Windows-x86.exe' C:\Miniconda.exe; echo "Done"}
  - cmd: C:\Miniconda.exe /S /D=C:\Miniconda
  - ps: ls C:\Miniconda/Scripts
  # We need to do this first as other commands may not work with older versions of conda.
  - C:\Miniconda\Scripts\conda update -yq conda
  - C:\Miniconda\Scripts\conda install -yq -c conda pytest requests conda-env
  - ps: if($env:PYTHON -eq "2.7"){C:\Miniconda\Scripts\conda install -yq ndg-httpsclient}
  # Other packages
  - C:\Miniconda\Scripts\conda install -n conda atlas numpy scipy matplotlib nose dateutil pandas statsmodels pip --quiet

install:


  - C:\Miniconda\python.exe setup.py install
  - C:\Miniconda\Scripts\conda info -a
  - C:\Miniconda\Scripts\conda list

test_script:
  - C:\Miniconda\Scripts\py.test -vvv
  - IF %PYTHON%==2.7 IF DEFINED PYCOSAT C:\Miniconda\Scripts\conda install -f menuinst ipython-notebook & echo "Finished menuinst test"
