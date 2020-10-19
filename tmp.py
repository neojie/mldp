import os
def load_paths(inputpath,level='vasp',output=False):
    """
    load paths and modify path to level required
    
    inputpath: input path, a file name storing input paths, str
    level: level could be vasp, recal, deepmd
    output: write output as a single file
    
    return
    ------
    paths: list of str
    
    Note
    -----
    path should be:
    xx/recal/deepmd*, xx/recal, xx
    or
    xx/recal/deepmd*  relax_step, relax_paths
    
    xx is vasp level
    deepmd* could be deepmd, deepmd_relax, deepmd_relax2, deepmd_relax, etc.
    

    """
    fp          = open(inputpath,'r')
    folders_org = fp.readlines()
    paths     = []
    relax_paths = []
    relax_step  = []
    fp.close()
    relax = False
    for i in range(len(folders_org)):
#        print(i)
        #skip #xxx and empty
        if folders_org[i].split() == [] or '#' == folders_org[i].split()[0] or folders_org[i].split()[0] == '\n':
            pass
        else:
            # remove \n in the end
            tmp = folders_org[i].replace('\n','')
            # remove # 
            tmp = tmp.split('#')[0]
            # remove " and ,                     
            tmp= tmp.replace('"','').replace(',','').split()
            
            if len(tmp) == 3:
                relax = True
                tmp2 = _load_paths_helper(tmp[0],level)
                if not(tmp2 in paths):
                    paths.append(tmp2)
                    
                if tmp[2] == 'self':
                    relax_paths.append(paths[-1])
                else:
                    relax_paths.append(tmp[2])
                relax_step.append(int(tmp[1]))
                
            elif len(tmp) == 1:
                tmp2 = _load_paths_helper(tmp[0],level)
                if not(tmp2 in paths):
                    paths.append(tmp2)
            else:
                print('????????',folders_org[0],'length is not 3 or 1, skip!')
                
    
    if relax_paths != []:
        paths = [paths,relax_step, relax_paths]
    
    if output:
        from datetime import date
        today = date.today()
        out = open('paths-{0}-{1}-{2}'.format(today.year,today.month,today.day),'w')
        if relax is True:
            for path in paths[0]:
                out.write(path)
                out.write('\n')
        else:
            for path in paths:
                out.write(path)
                out.write('\n')
        out.close()            
    return paths


def _load_paths_helper(path,level):
    """
    if '/recal' exist
    use recal to divide the path
    if recal does not exist
    must be vasp level
    """
    if '/recal' in path:
        vasp_path = path.split('/recal')[0]
    else:
        vasp_path = path
    
    if level == 'vasp':
        return vasp_path
    if level == 'recal':
        return os.path.join(vasp_path,'recal')
    if 'deepmd' in level:
        return os.path.join(vasp_path,level)
        

