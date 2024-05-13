import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="gzjwst",
    version="0.0.1",
    author="Galaxy Zoo collaboration",
    author_email="walmsleymk1@gmail.com",
    description="Pipeline for GZ JWST cutouts",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/mwalmsley/gz-jwst",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License (GPL)",
        "Operating System :: OS Independent",
        "Development Status :: 4 - Beta"
    ],
    packages=setuptools.find_packages(),
    python_requires=">=3.9",  # bumped to 3.9 for typing
    install_requires=[
        'matplotlib',
        'seaborn',
        'jupyterlab',
        'tqdm',
        'pillow',
        'numpy',
        'pandas',
        'scipy',
        'astropy', 
        'matplotlib',
        'pyarrow',  # to read parquet, which is very handy for big datasets
        'setuptools',  # no longer pinned
    ]
)
