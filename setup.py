import sys
from distutils.core import setup
import codecs
try:
    codecs.lookup('mbcs')
except LookupError:
    ascii = codecs.lookup('ascii')
    func = lambda name, enc = ascii: {True: enc}.get(name=='mbcs')
    codecs.register(func)

DESCRIPTION = """\
One-dimensional modeling for land surface subsidence in Mexico City
"""

# setup(name="subtim", version="0.1", description="One-dimensional modeling for land surface subsidence in Mexico City",
#     author="Emily Pease", packages=["subtim"])

def run():
    setup(name="subpydence",
          version="0.1",
          description="One-dimensional modeling for land surface subsidence in Mexico City",
          author="Emily Pease",
          packages=["subpydence"],
          author_email = 'emilypease@utexas.edu'
          )
if __name__ == "__main__":
    run()