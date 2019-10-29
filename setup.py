import setuptools

with open("README.md", "r") as fh:
	long_description = fh.read()

setuptools.setup(
	name="SplitP",
	version="0.0.3",
	author="Joshua Stevenson",
	author_email="joshua.stevenson@utas.edu.au",
	description="SplitP",
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
)