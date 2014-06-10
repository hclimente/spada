from distutils.core import setup, Extension

module1 = Extension('C_functions',
                    #include_dirs = ['/home/emre/lib/Python-2.4.4/Include/', '/home/emre/lib/Python-2.4.4/', '/home/emre/lib/mysql-5.0.51a/include/'],
                    include_dirs = ['/usr/include/python2.5/', '/home/emre/lib/mysql-5.0.51a/include/'],
                    libraries = ['mysqlclient'],
                    #library_dirs = [ '/home/emre/lib/mysql-5.0.51a/libmysql' ],
                    library_dirs = [ '/usr/lib64/mysql' ],
                    sources = ['C_functions.cpp'])

setup (name = 'C_functions', 
       version = '1.2', 
       description = 'Biana compiled functions', #'This is a C package for optimization of unification procedures',
       ext_modules = [module1])

