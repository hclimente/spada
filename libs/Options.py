#!/soft/devel/python-2.7/bin/python

import argparse

class Options(object):
    _instance = None
    def __new__(cls, *args, **kwargs):
        if not cls._instance:
            cls._instance = super(Options, cls).__new__(
                                cls, *args, **kwargs)
        return cls._instance


    def parse_arguments(*args, **kwds):
        parser = argparse.ArgumentParser(description = "kks", epilog= "kk")
        parser.add_argument('-f','--config-file', dest='config_file', action = 'store', default='Parameters.cfg',
                            help = 'Input file with the parameters')
        options=parser.parse_args()
        return options


if __name__ == '__main__':
    s1=Options()
    s2=Options()
    if(id(s1)==id(s2)):
        print "Same"
    else:
        print "Different"