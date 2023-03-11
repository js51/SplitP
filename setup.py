import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

def readfile(filename):
    with open(filename,  encoding='utf-8') as f:
        return f.read()

setuptools.setup(
    name="splitp",
    version=readfile("VERSION").strip(),
    author="Joshua Stevenson",
    author_email="joshua.stevenson@utas.edu.au",
    description="Phylogenetic tools and methods involving splits and matrix rank",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/js51/SplitP",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
    package_data={'splitp': ['*.txt']},
    include_package_data=True,
)
