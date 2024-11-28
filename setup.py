import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="HCSIM",
    version="1.0.0",
    author="Xikang Feng",
    author_email="fxk@nwpu.edu.cn",
    maintainer="Sisi Peng",
    maintainer_email="sisipeng@mail.nwpu.edu.cn",
    description="HCSIM: A Single-Cell Genomics Simulator with Haplotype-Specific Copy Number Annotation.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/xikanfeng2/HCSIM",
    project_urls={
        "Bug Tracker": "https://github.com/xikanfeng2/HCSIM/issues",
    },
    classifiers = [
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    keywords=[
        'cancer',
        'single-cell',
        'DNA',
        'copy-number',
        'haplotype-specific',
        'simulator',
    ],
    package_dir={"": "src"},
    packages=setuptools.find_packages(where="src"),
    python_requires=">=3.6",
    install_requires=[
        'pandas>=0.23.4',
        'numpy>=2.1.0',
        'matplotlib>=3.0.2',
        'networkx>=3.2.1',
    ],
    entry_points={
        'console_scripts': [
            'hcsim=hcsim.bin.hcsim_main:main',
        ],
    },
)