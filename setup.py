try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup
import re
v_file = "case_control_filter/version.py"
v_line = open(v_file, "rt").read()
v_re = r"^__version__ = ['\"]([^'\"]*)['\"]"
match = re.search(v_re, v_line, re.M)
if match:
    verstr = match.group(1)
else:
    raise RuntimeError("Unable to find version string in {}.".format(v_file))

setup(
    name="case_control_filter",
    packages=["case_control_filter"],
    version=verstr,
    description="TODO - Project description",
    author="David A. Parry",
    author_email="david.parry@ed.ac.uk",
    url='https://github.com/david-a-parry/case_control_filter',
    download_url='https://github.com/david-a-parry/case_control_filter/archive/{}.tar.gz'.format(verstr),
    license='MIT',
    install_requires=['pysam>=0.14'],
    scripts=["bin/case_control_filter"],
    include_package_data=True,
    classifiers=[
        "Programming Language :: Python :: 3",
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Developers",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        ],
)
