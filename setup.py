from distutils.core import setup
from setuptools import find_packages

setup(
        name='hamster-tracking',
        version='1.0.0',
        license='gpl-3.0',
        description = 'HAMSTER: a Heat And MoiSture Tracking framEwoRk for tracking air parcels through the atmosphere',
        author = 'Jessica Keune, Dominik Schumacher, Diego G. Miralles',
        author_email = 'jessica.keune@ugent.be',
        url = 'https://github.com/h-cel/hamster',
        download_url ='https://github.com/h-cel/hamster/archive/v1.0.tar.gz',
        keywords = ['Lagrangian models', 'moisture tracking', 'precipitation recycling', 'origins of precipitation',
                    'origins of heat', 'source regions of moisture', 'source regions of heat', 'air parcel tracking', 
                    'land--atmosphere interactions'],   # Keywords
        # packages=['class4gl'],
        install_package_data = True,
        packages=find_packages("."),
        install_requires=['numpy','pandas'],
        long_description=open('README.md').read(),
        classifiers=[
                'Development Status :: 4 - Beta',
                'Intended Audience :: Atmospheric sciences, Hydrology',
                'Topic :: Tracking the origins of heat and moisture in the atmosphere',
                'License :: gpl-3.0',   
                'Programming Language :: Python :: 3',      
                # 'Programming Language :: Python :: 3.6',
                # 'Programming Language :: Python :: 3.7',
              ],
)
