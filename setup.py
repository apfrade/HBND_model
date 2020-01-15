import setuptools

with open("README.md", "r") as readme_file:
    long_description = readme_file.read()

requirements = ["pandas >= 0.2", "joblib >= 0.1", "numpy >= 1"]

setuptools.setup(
    name='HBNDmodel',
    version="0.0.8",
    author='Andre Frade',
    author_email="andre.frade@hertford.ox.ac.uk",
    description='HBND predictive model package',
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/apfrade/HBNDmodel.git",
    packages=setuptools.find_packages(),
    install_requires= requirements,
    python_requires='>=3',
    package_data={
    'HBNDmodel': ["data/descriptor_names", "data/hbnd_model.sav", "data/scaler_model.bin" ],
    },

    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
],  
)


