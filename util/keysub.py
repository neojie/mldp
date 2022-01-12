import argparse
from subprocess import call
parser = argparse.ArgumentParser(description="manipulate INCAR or other files")
parser.add_argument("--file","-f", type=str,default='INCAR', help="| file to modify")
parser.add_argument("--key_sub","-ks",nargs= 2, help="| repalce key. \n e.g., 'A = 1', 'A', 'B'=> 'B = 1'")
parser.add_argument("--key_val_sub","-kvs",nargs= 2, help="| repalce value of key. e.g., 'A = 1', 'A', '2'=> 'A = 2'")
parser.add_argument("--key_line_sub","-kls",nargs= 2, help="| repalce whole line containing key, e.g., 'A = 1', 'A', 'B=1'=> 'B=1'")
parser.add_argument("--add_line","-al",nargs= 1, help="| add a new line")
args = parser.parse_args()


if args.key_val_sub:
    com =  '/{0}/s/\([0-9].\+\)/{1}/g'.format(args.key_sub[0],args.key_sub[1])
    com = "sed -i -e '{0} ' {1}".format(com,args.file)
    print(com)
    call(com,shell=True) # sed to replace numbers after SIGMA, 

# sed -i '/LWAVE/c\   LWAVE        = .FALSE.' INCAR # sed to replace the whole line    
if args.key_line_sub:
    com =  '/{0}/c\{1}'.format(args.key_line_sub[0],args.key_line_sub[1])
    com = "sed -i '{0} ' {1}".format(com,args.file)
    print(com)
    call(com,shell=True) # sed to replace numbers after SIGMA, 

if args.key_sub:
    print('TODO')

if args.add_line:
    com = 'cat {0} >> {1}'.format(args.add_line, args.file)
    print(com)
    call(com,shell=True) 