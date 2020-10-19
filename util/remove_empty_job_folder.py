import os
import shutil
files = os.listdir()
kept    = 0 
removed = 0
for file in files:
    if os.path.isdir(file): 
        if os.path.exists(os.path.join(file,'OUTCAR')):
           kept +=1 
        else:
           shutil.rmtree(file)
           removed +=1
    else:
        pass
print('{0} jobs have outcar'.format(kept))
print('{0} jobs removed'.format(removed))
