from setuptools import setup, find_packages
from pathlib import Path

this_directory = Path(__file__).parent
long_description = (this_directory / "README.md").read_text(encoding="utf-8")

setup(
    name="MDInteractions",
    version="1.0.0",
    packages=find_packages(),
    install_requires=[
        "MDAnalysis",
        "numpy",
        "pandas"
    ],
    entry_points={
        'console_scripts': [
            'protein_interactions = MDInteractions.cli:main',  
            'mean_distance = MDInteractions.mean_cli:main'  
        ]
    },
    author="Anjani Fowdar and Darren P Martin",
    author_email="anjanifowdar@gmail.com", 
    description="A python package for the distance-based analysis of intra- and inter-protein interactions",
    long_description=long_description,
    long_description_content_type="text/markdown", 
    url="https://github.com/afowdar/MDInteractions", 
    python_requires=">=3.10", 
    classifiers=[
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Programming Language :: Python :: 3.12",
        "Operating System :: OS Independent",
        "License :: OSI Approved :: MIT License",
        "Development Status :: 5 - Production/Stable",  
        "Intended Audience :: Science/Research",
        "Intended Audience :: Developers",
        "Topic :: Scientific/Engineering :: Chemistry",
        "Topic :: Scientific/Engineering :: Bio-informatics",
        "Topic :: Software Development :: Libraries :: Python Modules",
    ],
)
