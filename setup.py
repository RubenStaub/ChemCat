import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="ChemCat",
    version="1.0",
    author="Ruben Staub",
    author_email="ruben.staub@ens-lyon.fr",
    description="Concatenate chemical structures easily",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/RubenStaub/ChemCat",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU Lesser General Public License v3 (LGPLv3)",
        "Operating System :: OS Independent",
        "Topic :: Scientific/Engineering :: Chemistry",
        "Intended Audience :: Science/Research",
    ],
#    python_requires='>=3.5', # python2 support is left to the motivated user...
    install_requires=['ase', 'numpy']
)

