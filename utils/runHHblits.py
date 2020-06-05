# -*- coding:utf-8 -*-
import os

class HHblits():
    def runHHblits(self,fastapath,outpath):
        names = [name for name in os.listdir(fastapath) if os.path.isfile(os.path.join(fastapath + '//', name))]
        for each_item in names:
            pdb_id = each_item.split('.')[0]
            '''
            database used: uniprot20
            link: http://wwwuser.gwdg.de/~compbiol/data/hhsuite/databases/hhsuite_dbs/
            '''
            cmd = '/home/lznenu/hh-suite/build/bin/hhblits -i '+ fastapath + '/' + each_item + ' -ohhm ' + outpath + '/' + pdb_id + '.hhm -d /home/lznenu/hh-suite/databases/uniprot20_2016_02/uniprot20_2016_02 -v 0 -maxres 40000 -Z 0'
            os.system(cmd)
      
      
if __name__ == '__main__':

    '''
    samples:
    fastapath = '/home/liuzhe002/HHblits/HHFasta'
    outpath = '/home/RaidDisk/liuzhe002/HHResult'
    You can also check the structure and format of used fasta files in my folder : /home/liuzhe002/HHblits/HHFasta
    Warning : the permissions issue has not been resolved, please use the filepath under /home/RaidDisk/ as your outpath
    '''
    
    fastapath = '/home/lznenu/ss_topo/test'
    outpath = '/home/lznenu/ss_topo/ali_HHResults'
    
    hh = HHblits()
    hh.runHHblits(fastapath, outpath)
