import os
try:
    from setuptools import setup, find_packages
except ImportError:
    from distutils.core import setup, find_packages


def read(fname):
    # Utility function to read the README file.
    # Used for the long_description.  It's nice, because now 1) we have a top level
    # README file and 2) it's easier to type in the README file than to put a raw
    # string in below ...
    return open(os.path.join(os.path.dirname(__file__), fname)).read()


setup(
    name="pinfer",
    version="v0.7.2b",
    author="Nick Fyson",
    author_email="nick@fyson.net",
    description="Inference of ancestral Protein Interaction Networks (PINs).",
    long_description=read('README.md'),
    install_requires=read('requirements.txt'),
    license="BSD",
    keywords='PPI protein modelling interaction',
    url='https://github.com/nickfyson/pinfer',
    packages=find_packages(),
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: BSD License',
        'Natural Language :: English',
        'Programming Language :: Python',
        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
        'Operating System :: Microsoft :: Windows',
        'Operating System :: POSIX',
        'Operating System :: Unix',
        'Operating System :: MacOS',
    ],
)
