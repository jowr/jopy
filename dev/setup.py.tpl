import os
import sys

try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup


if sys.argv[-1] == 'publish':
    os.system('python setup.py sdist upload')
    sys.exit()

readme = open('README.rst').read()
doclink = """
Documentation
-------------

The full documentation is at http://jopy.rtfd.org."""
history = open('HISTORY.rst').read().replace('.. :changelog:', '')

setup(
    name='jopy',
    version='{{ version }}',
    description='My swiss army knife for plotting and Modelica simulations in Python',
    long_description=readme + '\n\n' + doclink + '\n\n' + history,
    author='{{ author }}',
    author_email='{{ email }}',
    url='https://github.com/jowr/jopy',
    packages=[
        'jopy',
    ],
    package_dir={'jopy': 'jopy'},
    include_package_data=True,
    install_requires=[{% for pkg in cus_pkgs %}
      '{{ pkg -}}',{% endfor %}{% for pkg in pip_cus_pkgs %}
      '{{ pkg -}}',{% endfor %}
    ],
    license='MIT',
    zip_safe=False,
    keywords='jopy',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: Implementation :: PyPy',
    ],
)
