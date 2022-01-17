import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name='pemfc_model_lufire',
    version="0.0.1",
    author="Lukas Feierabend",
    author_email="lukas.feierabend@gmail.com",
    description="PEM Fuel Cell Stack Simulation Model",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/ZBT-Tools/PEMFC-Model",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)