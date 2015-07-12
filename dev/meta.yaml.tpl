package:
  name: jopy
  version: {{ version }}

source:
  git_tag: {{ tag }}
  git_url: https://github.com/jowr/jopy.git

  # If this is a new build for the same version, increment the build
  # number. If you do not include this key, it defaults to 0.
  # number: 1

build:
  script: python setup.py install

requirements:
  build:
    - python
    - setuptools{% for pkg in dev_pkgs %}
    - {{ pkg -}}
{% endfor %}{% for pkg in pip_dev_pkgs %}
    - {{ pkg -}}
{% endfor %}

  run:
    - python{% for pkg in cus_pkgs %}
    - {{ pkg -}}
{% endfor %}{% for pkg in pip_cus_pkgs %}
    - {{ pkg -}}
{% endfor %}

test:

  imports:
    - jopy

  commands:
    # You can put test commands to be run here.  Use this to test that the
    # entry points work.
    nosetests --with-coverage --cover-package=jopy --logging-level=INFO

  # You can also put a file called run_test.py in the recipe that will be run
  # at test time.

  requires:{% for pkg in dev_pkgs %}
    - {{ pkg -}}
{% endfor %}{% for pkg in pip_dev_pkgs %}
    - {{ pkg -}}
{% endfor %}

about:
  home: http://www.jorrit.org
  license: MIT License
