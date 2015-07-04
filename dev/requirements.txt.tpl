package:
  name: jopy
  version: {{ version }}

source:
  git_tag: {{ tag }}
  git_url: https://github.com/jowr/jopy.git

  # If this is a new build for the same version, increment the build
  # number. If you do not include this key, it defaults to 0.
  # number: 1

requirements:
  build:
    - python
    - setuptools{% for pkg in pkgs %}
    - {{ pkg -}}
{% endfor %}

  run:
    - python{% for pkg in pkgs %}
    - {{ pkg -}}
{% endfor %}

#test:
  # Python imports
  #imports:
  #  - jopy

  # commands:
    # You can put test commands to be run here.  Use this to test that the
    # entry points work.


  # You can also put a file called run_test.py in the recipe that will be run
  # at test time.

  # requires:
    # Put any additional test requirements here.  For example
    # - nose

about:
  home: http://www.jorrit.org
  license: MIT License
