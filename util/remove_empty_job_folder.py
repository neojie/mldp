import os
#import shutil
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--folder","-f",type=str,help="1-3 means dirs 1,2,3")
args   = parser.parse_args()

folder_range =[int(i) for i in args.folder.split('-')]    
folders = list(range(folder_range[0],folder_range[1]))
folders.append(folder_range[-1])
for i in folders: 
    try:
        print("remove ",i)
        os.removedirs(i)    
    except:
        print(i,"is not empty")
    
#files = os.listdir()
#kept    = 0 
#removed = 0
#for file in files:
#    if os.path.isdir(file): 
#        if os.path.exists(os.path.join(file,'OUTCAR')):
#           kept +=1 
#        else:
#           shutil.rmtree(file)
#           removed +=1
#    else:
#        pass
#print('{0} jobs have outcar'.format(kept))
#print('{0} jobs removed'.format(removed))
