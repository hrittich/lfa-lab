#!/usr/bin/env python
#
# LFA Lab - Library to simplify local Fourier analysis.
# Copyright (C) 2018-2020  Hannah Rittich
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

from setuptools import setup
from setuptools import Command
from setuptools.command.install import install
from distutils.command.build import build
from subprocess import check_call, CalledProcessError
import os
import sys
import distutils.spawn

if 'CMAKE_EXECUTABLE' in os.environ:
    CMAKE_EXE = os.environ['CMAKE_EXECUTABLE']
else:
    CMAKE_EXE = distutils.spawn.find_executable('cmake')

print('Using cmake executable: {}'.format(CMAKE_EXE))
check_call([CMAKE_EXE, '--version'])

class install_cmake(Command):

    user_options = [
        ('install-dir=', 'd', "directory to install to")
    ]

    def initialize_options(self):
        self.install_lib = None

    def finalize_options(self):
        # get attibute install_lib from install command and store into
        # install_lib variable
        self.set_undefined_options('install',
                                   ('install_lib', 'install_lib'))

        self.install_lib = os.path.abspath(self.install_lib)

    def run(self):
        check_call([CMAKE_EXE,
            '-DPYTHON_INSTALL_DIR={}'.format(self.install_lib),
            ])
        check_call(['make', 'install'])

    def get_outputs(self):
        outputs = []
        # Store all files installed by CMake
        # CMake stores all installed files in install_manifest.txt
        with open('install_manifest.txt') as in_fp:
            for line in in_fp:
                outputs.append(line.strip('\n\r'))

        print('Outputs returned by CMake: {}'.format(outputs))
        return outputs

class build_cmake(Command):
    def initialize_options(self):
        pass

    def finalize_options(self):
        pass

    def run(self):
        # run cmake
        try:
            check_call([CMAKE_EXE, '-DPYTHON_EXECUTABLE='+sys.executable, '.'])
        except CalledProcessError:
            # if cmake fails, try to run with a clean cache
            os.remove('CMakeCache.txt')
            check_call([CMAKE_EXE, '-DPYTHON_EXECUTABLE='+sys.executable, '.'])

class custom_build(build):
    sub_commands = [ ('build_cmake', lambda self: True) ] \
        + build.sub_commands

class custom_install(install):
    sub_commands = [ ('install_cmake', lambda self: True) ] \
        + install.sub_commands


setup(
    name='lfa-lab',
    version='0.4.0',
    cmdclass=dict(
        build_cmake=build_cmake,
        install_cmake=install_cmake,
        build=custom_build,
        install=custom_install,
        ),
    install_requires=[
        'numpy',
        'six'
    ]
    )
