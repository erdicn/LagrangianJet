import re
import os
from pathlib import Path
from data_analysis import *

def getAllCAData(base_dir = Path('.'),  
                 exclude_folders = {'OutCstm', "OutTenThou", 'Out1', 'Out2', 'tmp'},
                 rho = 1,
                 file_dim_to_SI = 1e-18, get_from_dir = False,
                 file_name_regex = "mass*.dat"):
    all_data_list = np.array([])
    all_dirs = None
    
    if get_from_dir:
        all_dirs = Path(get_from_dir).rglob(file_name_regex)
    else: 
        all_dirs = base_dir.rglob(file_name_regex)
    
    for file_path in all_dirs:
        parent_folder = file_path.parent.name
        if parent_folder in exclude_folders:
            continue
        try:
            if os.path.getsize(file_path) > 0:
                data = np.loadtxt(file_path, delimiter=" ")
                if data.size > 1:
                    all_data_list = np.append(all_data_list, data)
    
        except Exception as e:
            print(f"Skipped {file_path} due to error: {e}")
    return massToRadius(all_data_list*file_dim_to_SI, rho=rho)*2


def getDataFromFolder(base_dir = Path('.'),  
                      exclude_folders = {'OutCstm', "OutTenThou", 'Out1', 'Out2', 'tmp'},
                      get_from_dir = False,
                      file_name_regex = "radius*.dat"):
    all_data_list = np.array([])
    all_dirs = None
    
    if get_from_dir:
        all_dirs = Path(get_from_dir).rglob(file_name_regex)
    else: 
        all_dirs = base_dir.rglob(file_name_regex)
    
    valid_files = [
        file_path for file_path in all_dirs 
        if file_path.parent.name not in exclude_folders
    ]
    
    def extract_file_number(file_path):
        match = re.search(r'\d+', file_path.name)
        return int(match.group()) if match else float('inf')
    
    valid_files.sort(key=extract_file_number)
    
    collected_data = []
    
    for file_path in valid_files:
        try:
            if os.path.getsize(file_path) > 0:
                data = np.loadtxt(file_path, delimiter=" ")
                if data.size > 1:
                    collected_data.append(data)
                    
        except Exception as e:
            print(f"Skipped {file_path} due to error: {e}")
            
    # Combine everything at the end
    if collected_data:
        all_data_list = np.concatenate(collected_data)
    else:
        all_data_list = np.array([])
        
    # for file_path in all_dirs:
    #     parent_folder = file_path.parent.name
    #     if parent_folder in exclude_folders:
    #         continue
    #     try:
    #         if os.path.getsize(file_path) > 0:
    #             data = np.loadtxt(file_path, delimiter=" ")
    #             if data.size > 1:
    #                 all_data_list = np.append(all_data_list, data)
    
    #     except Exception as e:
    #         print(f"Skipped {file_path} due to error: {e}")
    return all_data_list