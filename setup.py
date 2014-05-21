from distutils.core import setup, Extension
import os, sys, glob

__version__ = '0.0.4'

setup(name = 'ratty2',
    version = __version__,
    description = 'Interfaces to an RFI monitoring spectrometer',
    long_description = 'Provides interfaces to a ROACH-based RFI monitoring spectrometer.',
    license = 'GPL',
    author = 'Jason Manley',
    author_email = 'jason_manley at hotmail.com',
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Developers',
        'Operating System :: OS Independent',
        'License :: OSI Approved :: GNU General Public License (GPL)',
        'Topic :: Scientific/Engineering :: Astronomy',
        'Topic :: Software Development :: Libraries :: Python Modules',
        ],
    requires=['h5py', 'pylab','matplotlib','numpy','corr'],
    provides=['ratty2'],
    package_dir = {'ratty2':'src'},
    packages = ['ratty2'],
    scripts=glob.glob('scripts/*'),
    data_files=[('/etc/ratty2', glob.glob('config_files/*')),('/etc/ratty2/cal_files', glob.glob('cal_files/*')),('/etc/ratty2/boffiles', glob.glob('boffiles/*'))]
)

