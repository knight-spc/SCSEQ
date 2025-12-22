# 用于解压压缩包
import zipfile
import os

def unzip(final_filename,save_path,user_id=0,task_id=0):
    zip_path=final_filename

    if not os.path.isdir(save_path):
        os.mkdir(save_path)
    save_path = save_path

    file=zipfile.ZipFile(zip_path)
    print('开始解压...')
    file.extractall(save_path)
    print('解压结束。')
    file.close()
def add_numbers(a, b):
    return a + b