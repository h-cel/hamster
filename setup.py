from distutils.core import setup
from setuptools import find_packages

setup(
        name='hamster-tracking',
        version='1.0.0',
        license='gpl-3.0',
        description = 'HAMSTER: a Heat And MoiSture Tracking framEwoRk for tracking air parcels through the atmosphere',
        author = 'Jessica Keune, Dominik L. Schumacher, Diego G. Miralles',
        author_email = 'jessica.keune@ugent.be',
        url = 'https://github.com/h-cel/hamster',
        download_url ='https://github.com/h-cel/hamster/archive/v1.0.0.tar.gz',
        keywords = ['Lagrangian models', 'moisture tracking', 'precipitation recycling', 'origins of precipitation',
                    'origins of heat', 'source regions of moisture', 'source regions of heat', 'air parcel tracking', 
                    'land--atmosphere interactions'],   # Keywords
        install_package_data = True,
        packages=find_packages("."),
        install_requires=['numpy','pandas','gzip','hdf5','netcdf4','argparse','calendar','os','fnmatch','imp','math','random','re','struct','sys','time','timeit','warnings','datetime','functools','dateutil'],
        long_description=open('README.md').read(),
        classifiers=[
                'Development Status :: 4 - Beta',
                'Intended Audience :: Atmospheric sciences, Hydrology',
                'Topic :: Tracking the origins of heat and moisture in the atmosphere',
                'License :: gpl-3.0',   
                'Programming Language :: Python :: 3',      
              ],
        package_dir={"": "src"},
        packages=setuptools.find_packages(where="src"),
        python_requires=">=3.6",
)
