from setuptools import setup, find_packages
import os

# プロジェクトルートのパスを取得
here = os.path.abspath(os.path.dirname(__file__))

# README.mdを読み込み
with open(os.path.join(here, "README.md"), "r", encoding="utf-8") as fh:
    long_description = fh.read()

# requirements.txtを読み込み
with open(os.path.join(here, "requirements.txt"), "r", encoding="utf-8") as fh:
    requirements = [line.strip() for line in fh if line.strip() and not line.startswith("#")]

setup(
    name="segene-analyzer",
    version="1.0.0",
    author="norio shinkai",
    description="SEdb Region Analyzer - Analyze super-enhancer regions using SEdb data",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/hamamoto-lab/SEgene/SEgene_analyzer",
    packages=find_packages(),
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.11",
        "Programming Language :: Python :: 3.12",
        "Programming Language :: Python :: 3.13",
    ],
    python_requires=">=3.11",
    install_requires=requirements,
    entry_points={
        'console_scripts': [
            'sedb-analyzer=segene_analyzer.cli:main',  # パッケージ名を変更
        ],
    },
    include_package_data=True,
)