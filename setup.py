"""
MethSCAn2 Setup Configuration
"""

from setuptools import setup, find_packages

# Read long description
with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

# Core dependencies
install_requires = [
    "numpy>=1.20.0",
    "scipy>=1.7.0",
    "pandas>=1.3.0",
    "anndata>=0.8.0",
    "scikit-learn>=1.0.0",
    "matplotlib>=3.5.0",
    "seaborn>=0.11.0",
    "umap-learn>=0.5.0",
    "leidenalg>=0.8.0",
    "python-igraph>=0.9.0",
    "joblib>=1.1.0",
    "tqdm>=4.62.0",
]

setup(
    name="methscan2",
    version="0.1.0",
    author="MethSCAn2 Team",
    author_email="your.email@example.com",
    description="Single-cell DNA methylation analysis toolkit",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/yourusername/methscan2",
    packages=find_packages(),
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
    ],
    python_requires=">=3.8",
    install_requires=install_requires,
    extras_require={
        "dev": [
            "pytest>=7.0",
            "pytest-cov>=3.0",
            "black>=22.0",
            "flake8>=4.0",
        ],
    },
    include_package_data=True,
    zip_safe=False,
)
