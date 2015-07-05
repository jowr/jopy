import jinja2
import os
import requests
import json
from distutils.version import LooseVersion #, StrictVersion
import codecs

template_dir = os.path.dirname(__file__)
target_dir = os.path.abspath(os.path.join(os.path.dirname(__file__),'..'))

bas_pkgs = ["conda-env", "conda-build", "binstar", "binstar-build", "jinja2"]
cus_pkgs = ["pip", "jinja2", "pyaml", "numpy", "scipy", "matplotlib", "nose", "dateutil", "pandas", "statsmodels", "pytest", "cython"]
pip_pkgs = ["coveralls", "nose-cov", "codecov", "coolprop", "texttable"]

loader = jinja2.FileSystemLoader(template_dir)
environment = jinja2.Environment(loader=loader)

target = 'travis.yml'
template = environment.get_template(os.path.join(template_dir,target+'.tpl'))
f = codecs.open(os.path.join(target_dir,"."+target),mode='wb',encoding='utf-8')
f.write(template.render(bas_pkgs=" ".join(bas_pkgs), cus_pkgs=" ".join(cus_pkgs), pip_pkgs=" ".join(pip_pkgs)))
f.close()

tags = {}
#r = requests.get('https://api.github.com/repos/coolprop/coolprop/tags')
r = requests.get('https://api.github.com/repos/jowr/jopy/tags')
if(r.ok):
    item = json.loads(r.text or r.content)
    for com in item:
        tags[com['name']] = com['commit']['sha']

#tag = sorted(tags.keys())[-1]
tag = "v0.0.1"
if tag[0]=='v': version = tag[1:]
else: version = tag

target = 'appveyor.yml'
template = environment.get_template(os.path.join(template_dir,target+'.tpl'))
f = codecs.open(os.path.join(target_dir,target),mode='wb',encoding='utf-8')
f.write(template.render(version=version, bas_pkgs=" ".join(bas_pkgs), cus_pkgs=" ".join(cus_pkgs), pip_pkgs=" ".join(pip_pkgs)))
f.close()

#target = 'meta.yaml'
#template = environment.get_template(os.path.join(template_dir,target+'.tpl'))
#tags = {}
#r = requests.get('https://api.github.com/repos/coolprop/coolprop/tags')
#if(r.ok):
    #item = json.loads(r.text or r.content)
    #for com in item:
        #tags[com['name']] = com['commit']['sha']

##tag = sorted(tags.keys())[-1]

##for tag in sorted(tags.keys()):
##    print tag
##     r = requests.get('https://api.github.com/repos/coolprop/coolprop/git/tags/'+tags[tag])
##     if(r.ok):
##         items = json.loads(r.text or r.content)
##         print str(items)

##def cmp(x,y): return LooseVersion(x).__cmp__(y)
##tag = sorted(tags.keys(),cmp=cmp)[-1]
#tag = sorted(tags.keys())[-1]
##from pkg_resources import parse_version
##>>> parse_version('1.4') > parse_version('1.4-rc2')
#if tag[0]=='v': version = tag[1:]
#else: version = tag

#f = codecs.open(os.path.join(target_dir,target),mode='wb',encoding='utf-8')
##f = open(name,mode='w')
#f.write(template.render(version=version,tag=tag,pkgs=pkgs))
#f.close()