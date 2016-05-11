#!/usr/bin/env python
# -*- coding: utf-8 -*-


try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup


with open('README.md') as readme_file:
    readme = readme_file.read()

with open('requirements.txt') as readme_file:
    requirements = readme_file.readlines()

test_requirements = [
    # TODO: put package test requirements here
]

setup(
    name='pinfer',
    version='0.6rc',
    description="Inference of ancestral Protein Interaction Networks (PINs).",
    long_description=readme + '\n\n',
    author="Nick Fyson",
    author_email='mail@nickfyson.co.uk',
    url='https://github.com/nickfyson/pinfer',
    packages=[
        'pinfer',
    ],
    package_dir={'pinfer':
                 'pinfer'},
    include_package_data=True,
    install_requires=requirements,
    license="BSD",
    zip_safe=False,
    keywords='pinfer',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: BSD License',
        'Natural Language :: English',
        'Programming Language :: Python :: 3.4',
    ],
    test_suite='tests',
    tests_require=test_requirements
)
