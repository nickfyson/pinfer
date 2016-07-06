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
    version="0.7rc",
    author="Nick Fyson",
    author_email="nick@fyson.net",
    description="Inference of ancestral Protein Interaction Networks (PINs).",
    long_description=read('README.md'),
    install_requires=read('requirements.txt'),
    license="BSD",
    keywords='PPI protein modelling interaction',
    url='https://github.com/nickfyson/pinfer',
    packages=['pinfer', 'pinfer.infer', 'pinfer.itree', 'pinfer.visualise', 'tests'],
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: BSD License',
        'Natural Language :: English',
    ],
)
