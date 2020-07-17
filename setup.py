import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="Stretchme", # Replace with your own username
    version="0.2",
    author='Paweł Dąbrowski-Tumański',
    author_email='p.dabrowski-tumanski@uw.edu.pl',
    description="A package for analysis of stretching of proteins.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/ilbsm/protein_stretching",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)
