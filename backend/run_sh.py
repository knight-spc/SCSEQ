import os
import subprocess
import time
from database import add_new_task
from database import update_item_sh_id



def check_slurm_job(job_id):
    """检查SLURM任务状态"""
    try:
        result = subprocess.run(['squeue', '-j', job_id], capture_output=True, text=True)
        return len(result.stdout.strip().split('\n')) <= 1
    except subprocess.CalledProcessError:
        return True  
    

def runsh(command, sh_path):

    sh_content = f"""
"""

    sh_file_path = os.path.join("backend/sh",sh_path)
    with open(sh_file_path, 'w') as sh_file:
        sh_file.write(sh_content)

    os.chmod(sh_file_path, 0o755)  

    try:
        result = subprocess.run(['sbatch', sh_file_path], check=True, capture_output=True, text=True)
        job_id = result.stdout.strip().split()[-1]
        print(f"任务已提交，ID: {job_id}")
        
        max_attempts = 720  
        attempt = 0
        while attempt < max_attempts:
            print("开始等待")
            if check_slurm_job(job_id):
                print("任务已完成")
                break
            time.sleep(10)  
            attempt += 1
            print(f"等待任务完成，尝试次数：{attempt}")
        
        if attempt >= max_attempts:
            raise Exception("任务执行超时")
        return job_id
    except subprocess.CalledProcessError as e:
        print(f"Error occurred: {e}")


def runsh_item(command, sh_path, param):

    sh_content = f"""
"""
    print(sh_content)
    sh_file_path = os.path.join("backend/sh",sh_path)
    with open(sh_file_path, 'w') as sh_file:
        sh_file.write(sh_content)

    os.chmod(sh_file_path, 0o755)  


    try:
        result = subprocess.run(['sbatch', sh_file_path], check=True, capture_output=True, text=True)
        job_id = result.stdout.strip().split()[-1]
        print(f"任务已提交，ID: {job_id}")
        print(           param["item_id"],
            param["task_type"],
            param["result_path"],
            job_id)
        print(param["parameter"])

        task_id = add_new_task(
            item_id=param["item_id"],
            task_type=param["task_type"],
            parameter=param["parameter"],
            result_path=param["result_path"],
            sh_id=job_id
        )

        if task_id:
            print(f"成功创建任务，任务ID为：{task_id}")
        else:
            print("创建任务失败")
        
        max_attempts = 720  
        attempt = 0
        while attempt < max_attempts:
            print("开始等待")
            if check_slurm_job(job_id):
                print("任务已完成")
                break
            time.sleep(10)
            attempt += 1
            print(f"等待任务完成，尝试次数：{attempt}")
        
        if attempt >= max_attempts:
            raise Exception("任务执行超时")
        return job_id,task_id
    except subprocess.CalledProcessError as e:
        print(f"Error occurred: {e}")


def runsh_item_update(command, sh_path, job_id):
    
    sh_content = f"""
"""
    print(sh_content)
    sh_file_path = os.path.join("backend/sh",sh_path)
    with open(sh_file_path, 'w') as sh_file:
        sh_file.write(sh_content)

    os.chmod(sh_file_path, 0o755)

    try:
        result = subprocess.run(['sbatch', sh_file_path], check=True, capture_output=True, text=True)
        job_id_new = result.stdout.strip().split()[-1]
        print(f"任务已提交，ID: {job_id_new}")

        if update_item_sh_id(job_id, job_id_new):
            print(f"成功更新项目 {job_id} 的sh_id为 {job_id_new}")
        else:
            print(f"更新项目 {job_id} 的sh_id失败")
        
        max_attempts = 720  
        attempt = 0
        while attempt < max_attempts:
            print("开始等待")
            if check_slurm_job(job_id_new):
                print("任务已完成")
                break
            time.sleep(10)
            attempt += 1
            print(f"等待任务完成，尝试次数：{attempt}")
        
        if attempt >= max_attempts:
            raise Exception("任务执行超时")
        return job_id_new
    except subprocess.CalledProcessError as e:
        print(f"Error occurred: {e}")
        return None

