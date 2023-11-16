import setuptools

with open("README.md", "r", encoding="utf-8") as fhand:
    long_description = fhand.read()

setuptools.setup(
    name="TeloVision",
    version="v0.3.2",
    author="Tim Verschuren",
    author_email="t.verschuren@wi.knaw.nl",
    description="TeloVision is a python package that determines the presence \
        of telomeres at the ends of scaffolds and generates a visualisation.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/TimVerschuren/TeloVision",
    project_urls={
        "Bug Tracker": "https://github.com/TimVerschuren/TeloVision/issues",
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    install_requires=["pandas", "Bio", "plotly", "python-levenshtein"],
    packages=setuptools.find_packages(),
    python_requires=">=3.10",
    entry_points={
        "console_scripts": [
            "telovision = TeloVision_src.cli:main",
        ]
    }
)
