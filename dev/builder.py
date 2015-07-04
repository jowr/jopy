import jinja2
import os
import requests
import json
from distutils.version import LooseVersion #, StrictVersion
import codecs

template_dir = os.path.dirname(__file__)
target_dir = os.path.abspath(os.path.join(os.path.dirname(__file__),'..'))

pkgs = ["numpy","matplotlib","scipy"]
#wheel>=0.22
coolprop
texttable

loader = jinja2.FileSystemLoader(template_dir)
environment = jinja2.Environment(loader=loader)
#Now, if you'd like to include /path/to/templates/includes/sidebar.html in the /path/to/templates/index.html template, you'd write the following in your index.html:
#{% include 'includes/sidebar.html' %}

target = 'meta.yaml'
template = environment.get_template(os.path.join(template_dir,target+'.tpl'))
tags = {}
r = requests.get('https://api.github.com/repos/coolprop/coolprop/tags')
if(r.ok):
    item = json.loads(r.text or r.content)
    for com in item:
        tags[com['name']] = com['commit']['sha']

#tag = sorted(tags.keys())[-1]

#for tag in sorted(tags.keys()):
#    print tag
#     r = requests.get('https://api.github.com/repos/coolprop/coolprop/git/tags/'+tags[tag])
#     if(r.ok):
#         items = json.loads(r.text or r.content)
#         print str(items)

#def cmp(x,y): return LooseVersion(x).__cmp__(y)
#tag = sorted(tags.keys(),cmp=cmp)[-1]
tag = sorted(tags.keys())[-1]
#from pkg_resources import parse_version
#>>> parse_version('1.4') > parse_version('1.4-rc2')
if tag[0]=='v': version = tag[1:]
else: version = tag



f = codecs.open(os.path.join(target_dir,target),mode='wb',encoding='utf-8')
#f = open(name,mode='w')
f.write(template.render(version=version,tag=tag,pkgs=pkgs))
f.close()