__author__ = 'spark'

import os
import subprocess

__all__ = [
    "VariantCaller",
]

class VariantCaller():
    def __init__(self, home_dir, refgenome_filepath=None, inputbam_filepath=None, data_outdir='data/out/', tmpdir='data/results/tmp'):
        if (refgenome_filepath != None):
            self._refgenome_filepath = refgenome_filepath
        if (inputbam_filepath != None):
            self._inputbam_filepath = inputbam_filepath
        self._data_outdir = data_outdir
        self._tmpdir = tmpdir
        self._home_dir = home_dir

    def callVariants(self):
        try:
            self._callVariants(self._refgenome_filepath, self._inputbam_filepath, 'raw_variants.vcf')
        except Exception as e:
            print(e)

    def _callVariants(self, refgenome_filepath, inputbam_filepath, outputvcf_filename):
        try:
            subprocess.call(['java', '-Xmx4g', '-Djava.io.tmpdir='+self._tmpdir ,'-jar', self._home_dir+'toolset/GATK/gatk-package-4.6.0.0-local.jar',
                             'HaplotypeCaller',
                             '-R', refgenome_filepath,
                             '-I', inputbam_filepath,
                             '-ploidy', '1',
                             '-O', self._data_outdir + outputvcf_filename])
        except Exception as e:
            print(e)
