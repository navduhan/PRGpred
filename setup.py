import setuptools 
from pathlib import Path
this_directory = Path(__file__).parent
long_description = (this_directory / "README.md").read_text()
setuptools.setup(
    name='PRGpred',
    version='0.1',
    description ="PRGpred: Plant resistance genes prediction",
    url="https://bioinfo.usu.edu/PRGpred",
    author="Naveen Duhan",
    author_email = "naveen.duhan@usu.edu",
    license='MIT',
    packages = setuptools.find_packages(),
    include_package_data=True,
    long_description=long_description,
    long_description_content_type='text/markdown',
    entry_points={
            'console_scripts': [
                    'PRGpred = PRGpred.__main__:main',
            ]
    },
#     install_requires=[line.rstrip() for line in open("requirements.txt", "rt")],
install_requires = ['numpy','pandas','keras','sklearn','xlrd', 'openpyxl','h5py', 'biopython','scikit-learn','scipy','tensorflow'],

 classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: Unix",
        
        # Indicate who your project is intended for
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',


        # Specify the Python versions you support here. In particular, ensure
        # that you indicate whether you support Python 2, Python 3 or both.
        
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9'
    ],
    package_data={'': ['data/*.h5']},
    python_requires='>=3.8'
)