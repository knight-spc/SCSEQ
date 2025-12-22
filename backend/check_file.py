import os

def check_files(path):
    truefiles = {'barcodes.tsv.gz', 'features.tsv.gz', 'matrix.mtx.gz'}
    dirs = []
    for item in os.scandir(path):
        if item.is_file():
          print(item.name)
          if item.name != ".DS_Store":
            return False
        elif item.is_dir():
          dirs.append(item.path)

    for dir in dirs:
        subdir = ''
        chifiles = []
        for item in os.scandir(dir):
            if item.name != 'filtered_feature_bc_matrix':
                if item.name != ".DS_Store":
                    return False
            else:
                subdir = item.path
                
        for item in os.scandir(subdir):
            if item.name != ".DS_Store":
                chifiles.append(item.name)
        chifiles = set(chifiles)
        if chifiles != truefiles:
            print(chifiles)
            return False
    print("ok老铁没毛病！")
    return True



def copy_file_with_number(file_path, number):
    """
    复制文件并在文件名后添加数字
    Args:
        file_path (str): 原文件的完整路径
        number (int): 要添加的数字
    Returns:
        str: 新文件的路径
    """
    import os
    import shutil
    
    # 获取文件目录和文件名
    directory = os.path.dirname(file_path)
    filename = os.path.basename(file_path)
    
    # 分离文件名和扩展名
    name, ext = os.path.splitext(filename)
    
    # 创建新文件名
    new_filename = f"{name}_{number}{ext}"
    new_filepath = os.path.join(directory, new_filename)
    
    # 复制文件
    try:
        shutil.copy2(file_path, new_filepath)
        return new_filepath
    except Exception as e:
        print(f"复制文件失败: {e}")
        return None