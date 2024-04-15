import sys
from distutils.core import setup

if sys.version_info < (3, 1):
    sys.exit('Sorry, DIPAN requires Python >= 3.1 (Tested with Python 3.6.15)')

setup(name = 'DIPAN',
      version = "1.0", 
      description='DIPAN searches for neoantigens derived from intronic polyadenylation events found in tumor transcriptomes.',
      url='https://github.com/YY-TMU/DIPAN',
      scripts=['DIPAN.sh',
               'script/IPAFinder_DetectIPA.py',
               'script/build_transcript.py',
               'script/extract_sequence.py',
               'script/filter_normal_proteome.py'],
      install_requires=['HTseq == 0.9.1',
                        'numpy == 1.19.5',
                        'pandas == 1.1.5',
                        'pyfasta == 0.5.2',
                        'scipy == 1.5.4'],
      python_requires='~=3.6',
      )