#!/usr/bin/env python

###############################################################################
#                                                                             #
#    This program is free software: you can redistribute it and/or modify     #
#    it under the terms of the GNU General Public License as published by     #
#    the Free Software Foundation, either version 3 of the License, or        #
#    (at your option) any later version.                                      #
#                                                                             #
#    This program is distributed in the hope that it will be useful,          #
#    but WITHOUT ANY WARRANTY; without even the implied warranty of           #
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            #
#    GNU General Public License for more details.                             #
#                                                                             #
#    You should have received a copy of the GNU General Public License        #
#    along with this program. If not, see <http://www.gnu.org/licenses/>.     #
#                                                                             #
###############################################################################


import os
from setuptools import setup, find_packages


def version():
    setupDir = os.path.dirname(os.path.realpath(__file__))
    versionFile = open(os.path.join(setupDir, 'phyloutils', 'VERSION'), "r")
    _version = versionFile.readline().strip()
    versionFile.close()
    return _version


if __name__ == '__main__':

    dirName = os.path.dirname(__file__)
    if dirName and os.getcwd() != dirName:
        os.chdir(dirName)

    # packages=find_packages(),
    setup(
        name='phyloutils',
        version=version(),
        author='Shengwei Hou',
        author_email='housw2010@gmail.com',
        maintainer='Shengwei Hou',
        maintainer_email='housw2010@gmail.com',
        packages=['phyloutils', 'phyloutils.subcommands'],
        scripts=['bin/phyloutils'],
        url='http://pypi.python.org/pypi/phyloutils/',
        license='GPLv3',
        description='A set of scripts for phylogenetic manipulations.',
        classifiers=['Development Status :: 4 - Beta',
                     'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
                     'Natural Language :: English',
                     'Programming Language :: Python :: 3.7',
                     'Topic :: Scientific/Engineering :: Bio-Informatics',],
        install_requires=[
                          "biopython >= 1.68",
                          "ftputil >= 3.4",
        ],
        include_package_data=True,
        zip_safe=False,
        test_suite='nose.collector',
        tests_require=['nose'],
    )
