import os

__version__ = "1.0.0"


pkgname = "telofinder"

try:
    import easydev

    telofinder_path = easydev.get_package_location(pkgname)
    sharedir = os.sep.join([telofinder_path, pkgname, "data"])

    def get_data(filename):
        filename = sharedir + os.sep + filename
        if os.path.exists(filename):
            return filename
        else:
            raise IOError("unknown file {}".format(filename))


except:
    print('please install easydev with the command "pip install easydev"')
