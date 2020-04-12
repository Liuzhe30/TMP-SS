# -*- coding:utf-8 -*-
'''
	File Name：     runHHblits
	Description :   generate HHblits features of fasta files
	Author :        Liu Zhe
	date：          2020/2/21
'''
import os

class HHblits():
    def runHHblits(self,fastapath,outpath):
        names = [name for name in os.listdir(fastapath) if os.path.isfile(os.path.join(fastapath + '//', name))]
        for each_item in names:
            pdb_id = each_item.split('.')[0]
            '''
            database used: pdb70 released on 2020/2/19 downloaded on 2020/2/20 
            link: http://wwwuser.gwdg.de/~compbiol/data/hhsuite/databases/hhsuite_dbs/pdb70_from_mmcif_200219.tar.gz 
            '''
            cmd = '/home/ThirdPartTools/hh-suite/build/bin/hhblits -i '+ fastapath + '/' + each_item + ' -ohhm ' + outpath + '/' + pdb_id + '.hhm -d /home/RaidDisk/hh-suite/databases/pdb70_from_mmcif_200219/pdb70 -v 0 -maxres 40000 -Z 0'
            os.system(cmd)
      
      
if __name__ == '__main__':

    '''
    samples:
    fastapath = '/home/liuzhe002/HHblits/HHFasta'
    outpath = '/home/RaidDisk/liuzhe002/HHResult'
    You can also check the structure and format of used fasta files in my folder : /home/liuzhe002/HHblits/HHFasta
    Warning : the permissions issue has not been resolved, please use the filepath under /home/RaidDisk/ as your outpath
    '''
    
    fastapath = '/home/liuzhe002/2_gamma_turn/HHFasta'
    outpath = '/home/RaidDisk/liuzhe002/2_gamma_turn/HHResults'
    
    hh = HHblits()
    hh.runHHblits(fastapath, outpath)