import sys

__version__ = '$VERSION'

def get_version_message():
    # if I keep '$ VERSION' (without space) as is
    # the setup.py will replace it by the value set in setup
    # so the test become True even if craw is installed using setup.py
    if __version__ == '$' + 'VERSION':
        version = "NOT packaged, it should be a development version"
    else:
        version = __version__
    version_text = "craw {0} | Python {1}.{2}".format(version, sys.version_info.major, sys.version_info.minor)
    return version_text