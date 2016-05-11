import os
from setuptools import setup


def read(fname):
    # Utility function to read the README file.
    # Used for the long_description.  It's nice, because now 1) we have a top level
    # README file and 2) it's easier to type in the README file than to put a raw
    # string in below ...
    return open(os.path.join(os.path.dirname(__file__), fname)).read()


setup(
    name="pinfer",
    version="0.6",
    author="Nick Fyson",
    author_email="mail@nickfyson.co.uk",
    description="Inference of ancestral Protein Interaction Networks (PINs).",
    long_description=read('README.md'),
    install_requires=read('requirements.txt'),
    license="BSD",
    keywords='PPI protein modelling interaction',
    url='https://github.com/nickfyson/pinfer',
    packages=['pinfer', 'pinfer.infer', 'pinfer.itree', 'pinfer.visualise', 'tests'],
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: BSD License',
        'Natural Language :: English',
        'Programming Language :: Python :: 3.4',
    ],
)
