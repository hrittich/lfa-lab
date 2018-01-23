import numpy.distutils
import distutils.sysconfig as sc
import os

#from pprint import pprint
#pprint(sc.get_config_vars())

def cmake_list(l):
    return ';'.join(l)

def cmake_set_var(name, value, ctype, doc):
    print('set({0} "{1}" CACHE {2} "{3}")'.format(name, value, ctype, doc))

# print all NumPy include directories
ds = numpy.distutils.misc_util.get_numpy_include_dirs()
cmake_set_var('NUMPY_INCLUDES', cmake_list(ds), 'STRING', 
              'NumPy include directores')

cmake_set_var('NUMPY_FOUND', 'ON', 'BOOL', 'NumPy found')



cmake_set_var('PYTHON_INCLUDE_DIRS', sc.get_python_inc(), 'STRING',
              'Python includes')

python_lib = os.path.join(sc.get_config_var('LIBDIR'),
        sc.get_config_var('LDLIBRARY'))

cmake_set_var('PYTHON_LIB', python_lib, 'FILEPATH', 'Python libarary')

cmake_set_var('PYTHON_VERSION', sc.get_config_var('VERSION'), 'STRING',
              'Python Version')

cmake_set_var('PYTHON_DEFAULT_INSTALL_DIR', sc.get_python_lib(plat_specific=True),
    'PATH', 'Default Python directory which contains third party packages.')

