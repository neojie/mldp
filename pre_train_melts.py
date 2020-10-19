from shared_functions import load_paths
import os
import argparse
import numpy as np
parser = argparse.ArgumentParser()
### MUST SET###
parser.add_argument("--inputpath","-ip",help="input path file")

args = parser.parse_args()
cwd  = os.getcwd()


if args.inputpath:
    print("Check files in {0}  ".format(args.inputpath))
    inputpath = args.inputpath
    paths = load_paths(inputpath,'deepmd')


else:
    print("No folders point are provided. Use default value folders")
#    inputpath = os.path.join(cwd)
    paths = [cwd]
melts= []
for path in paths:
    set0 = os.path.join(path,'set.000')
    box = np.load(os.path.join(set0,'box.npy'))[0]
    if box[0] == box[4] and box[0] ==  box[-1]: #[12.805792,  0.      ,  0.      ,  0.      , 12.805792,  0.      ,0.      ,  0.      , 12.805792]
        melts.append(path)


with open('melts_to_train','w') as ft:
    ft.write('        "systems":      '+'[')
    for path in melts:
        ft.write('"'+path+'",')
        ft.write("\n")
    ft.write('],')
        
with open('melts_to_train2','w') as ft:
    for path in melts:
        ft.write(path)
        ft.write("\n")

