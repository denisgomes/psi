from setuptools import setup, find_packages


setup(
    name="psi",

    version="0.0.1",

    description=("Pipe Stress Infinity (PSI) - The pipe stress analysis and "
                 "design software."),

    # display on pypi
    long_description=open("README.rst", "r").read(),

    url="https://www.github.com/denisgomes/psi",

    author="Denis Gomes",

    author_email="denis.gomes@consultant.com",

    license="GPLv3",

    # advertise program attributes
    classifiers=[
        "Intended Audience :: End Users/Desktop",
        "Intended Audience :: Science/Research",
        "Intended Audience :: Developers",
        "License :: OSI Approved",
        "License :: OSI Approved :: GPL License",
        "Topic :: Scientific/Engineering",
        "Topic :: Scientific/Engineering :: Physics",
        "Topic :: Scientific/Engineering :: Visualization",
        "Topic :: Scientific/Engineering :: Artificial Intelligence",
        "Programming Language :: Python",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.5",
        ],

    keywords="piping pipe stress supports design engineering analysis",

    # excluded in build distributions, applies to packages only
    packages=find_packages(exclude=[
        "psi.docs", "psi.examples", "psi.tests"
        ]),

    # install from pypi, requirements.txt is for developers only
    install_requires=["scipy", "pint", "tqdm", "jinja2 >= 2.0"],

    package_data={
        "psi.data": ["pipes.csv", "materials.csv", "fluids.csv",
                     "insulation.csv"],
        "psi.data.spring": ["anvil_variable.csv", "anvil_constant.csv"],
        "psi.data.units": ["base.csv", "english.csv", "code_english.csv",
                           "si.csv"],
        },

    # MANIFEST.in works for source distributions only
    data_files=[("", ["LICENSE", "README.rst"])],

    # scripts= ,

    # tests
    test_suite="tests",

    entry_points={
        "console_scripts": [
            "psi = psi.bin.cli:main",
            ]
        }
    )
