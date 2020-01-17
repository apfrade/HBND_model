import setuptools

with open("README.md", "r") as readme_file:
    long_description = readme_file.read()

requirements = ["pandas >= 0.2", "joblib >= 0.1", "numpy >= 1"]

setuptools.setup(
    name='HBNDmodel',
    version="0.0.15",
    author='Andre Frade',
    author_email="andre.frade@hertford.ox.ac.uk",
    description='HBND predictive model package',
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/apfrade/HBNDmodel.git",
    packages=['HBNDmodel'],
    package_dir={'HBNDmodel': 'HBNDmodel'},
    package_data={'HBNDmodel':['HBNDmodel/data/descriptor_names', 'HBNDmodel/data/hbnd_model.sav', 'HBNDmodel/data/scaler_model.bin']},
    install_requires= requirements,
    python_requires='>=3',
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
    ],
    data_files=[(,[])]    
    include_package_data=True)


