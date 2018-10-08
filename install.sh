#!/bin/bash

python setup.py build
python setup.py install
python setup.py test
phyloutils -h
