import os
from datetime import datetime
import subprocess
import sqlite3
import pymysql

from flask import Flask, request, jsonify, send_file
import json

from werkzeug.utils import secure_filename
import time  
import socket
import shutil


import sys
import scanpy as sc
import celltypist
from celltypist import models
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

import pandas as pd
import numpy as np
from pandas.core.frame import DataFrame

from AICelltype import aicelltype

from unzip_file import unzip
from check_file import check_files, copy_file_with_number
from get_figure import pdf_to_png,pdf_to_svg
from read_csv import read_csv
from run_sh import runsh, runsh_item, runsh_item_update, check_slurm_job
import database

app = Flask(__name__)
app.config['UPLOAD_FOLDER'] = ''
app.config['DATA_FOLDER'] = ''
app.config['FIGURE_FOLDER'] = r''

app.config['MAX_CONTENT_LENGTH'] = 16 * 1024 * 1024 * 1024
host = ''
os.environ['QWEN_KEY'] = 'sk-'  # Add QWEN_KEY to environ, keep secret to yourself.



sqlite_save_path = ''
paths = ''

paths = ""
pathf = ''


def clear_directory(path):
    if not os.path.exists(path):
        print(f"路径 {path} 不存在")
        return
    
    for item in os.listdir(path):
        item_path = os.path.join(path, item)
        try:
            if os.path.isfile(item_path):
                os.remove(item_path)
            elif os.path.isdir(item_path):
                shutil.rmtree(item_path)
        except Exception as e:
            print(f"删除 {item_path} 时出错: {str(e)}")
    
    print(f"已清空目录 {path} 的内容")


@app.route('/api/create_item', methods=['POST'])
def create_item():
    print(request.headers)
    user_id = request.headers.get('user_id')
    name = request.form.get('name')       
    species = request.form.get('species')
    desc = request.form.get('desc')
    save_path = request.form.get('save_path')
    figure_path = request.form.get('figure_path')
    print("id:",user_id)
    print("输出：",name,species,desc)

    work_path = figure_path
    notes = desc

    item_id = database.add_item(user_id, name, save_path, work_path, species, notes)

    print("输入参数：",save_path,'\n',figure_path,'\n',species,'\n')
    Rpath = ""

    command = f'''Rscript {Rpath} -s "{save_path}" -f "{figure_path}" -p "{species}"'''
    print(command)
    sh_path = 'QC.sh'
    job_id = runsh(command, sh_path)
    print(f"任务{job_id}运行结束")

    print("create over")

    return jsonify(msg='成功', item_id=item_id), 200

@app.route('/api/get_user_items', methods=['GET'])
def get_user_items():
    try:
        user_id = request.headers.get('User-Id')
        if not user_id:
            return jsonify({
                "status": "error",
                "message": "未登录或登录已过期"
            }), 401

        try:
            conn = database.connect_db()
        except Exception as db_error:
            return jsonify({
                "status": "error",
                "message": f"数据库连接失败: {str(db_error)}",
                "error_type": "database_connection"
            }), 503 
        try:
            with conn.cursor() as cursor:
                sql = """SELECT item_id, item_name, species, notes, work_path, created_time 
                        FROM items WHERE user_id = %s 
                        ORDER BY item_id DESC"""
                cursor.execute(sql, (user_id,))
                items = cursor.fetchall()
                
                result = []
                for item in items:
                    result.append({
                        'id': item[0],
                        'name': item[1],
                        'species': item[2],
                        'notes': item[3],
                        'work_path': item[4],
                        'created_time': item[5].strftime('%Y-%m-%d %H:%M:%S') if item[5] else ''
                    })
                
                return jsonify({
                    "status": "success",
                    "items": result
                }), 200
        finally:
            conn.close()

    except Exception as e:
        print(f"获取用户项目失败: {str(e)}")
        return jsonify({
            "status": "error",
            "message": str(e)
        }), 500

@app.route('/api/get_item_tasks', methods=['GET'])
def get_item_tasks():
    item_id = request.headers.get('Item-Id')
    print("item_id:",item_id)
    if not item_id:
        return jsonify({"status": "error", "message": "缺少item_id"}), 400
    try:
        conn = database.connect_db()
        with conn.cursor() as cursor:
            cursor.execute("SELECT item_name FROM items WHERE item_id = %s", (item_id,))
            item_row = cursor.fetchone()
            project_name = item_row[0] if item_row else "未知项目"
        with conn.cursor() as cursor:
            sql = """SELECT task_id, type, parameter, time, sh_id FROM tasks WHERE item_id = %s ORDER BY time DESC"""
            cursor.execute(sql, (item_id,))
            tasks = cursor.fetchall()
        conn.close()
        tasks_copy = tasks

        cluster_map = {}
        for t in tasks:
            task_id, type_, parameter, time_, sh_id = t
            annotation_model = "未知"
            param_dict = {}
            if isinstance(parameter, str):
                try:
                    param_dict = json.loads(parameter.replace("'", '"'))
                except Exception:
                    import re
                    match = re.search(r'annotation_model[": ]+([a-zA-Z0-9_]+)', parameter)
                    if match:
                        annotation_model = match.group(1)
                else:
                    annotation_model = param_dict.get("annotation_model", "未知")
            elif isinstance(parameter, dict):
                annotation_model = parameter.get("annotation_model", "未知")
                param_dict = parameter

            if sh_id:
                is_finished = check_slurm_job(str(sh_id))
                status = "已完成" if is_finished else "运行中"
            else:
                status = "已完成"

            task_node = {
                "name": type_,
                "children": [
                    {"name": f"参数: {parameter}"},
                    {"name": f"提交时间: {time_.strftime('%Y-%m-%d %H:%M:%S') if time_ else ''}"},
                    {"name": f"状态: {status}"},
                ]
            }
            if annotation_model not in cluster_map:
                cluster_map[annotation_model] = []
            cluster_map[annotation_model].append(task_node)

        children = []
        for model, tasks in cluster_map.items():
            children.append({
                "name": model,
                "children": tasks
            })

        tree_data = {
            "name": project_name,
            "children": children
        }

        print("正在准备tree_data……")
        result = []
        for t in tasks_copy:
            task_id, type_, parameter, time_, sh_id = t

            if isinstance(parameter, str):
                parameter = parameter.strip('{}')
            elif isinstance(parameter, dict):
                parameter = ', '.join(f'{k}: {v}' for k, v in parameter.items())
            if sh_id:
                is_finished = check_slurm_job(str(sh_id))
                status = "已完成" if is_finished else "运行中"
            else:
                status = "已完成"  
            result.append({
                "task_id": task_id,
                "type": type_,
                "parameter": parameter,
                "time": time_.strftime('%Y-%m-%d %H:%M:%S') if time_ else "",
                "status": status
            })
        print()
        return jsonify({"status": "success", "tree": tree_data, "tasks": result}), 200
    except Exception as e:
        print(e)
        return jsonify({"status": "error", "message": str(e)}), 500       

@app.route('/api/delete_item', methods=['POST'])
def delete_item():
    try:
        item_id = request.headers.get('Item-Id')
        if not item_id:
            return jsonify({"status": "error", "message": "未登录或登录已过期"}), 401
        conn = database.connect_db()
        try:
            with conn.cursor() as cursor:
                sql_get_paths = "SELECT save_path, work_path FROM items WHERE item_id = %s"
                cursor.execute(sql_get_paths, (item_id,))
                paths = cursor.fetchone()
                
                if paths:
                    save_path, work_path = paths
                    for path in [save_path, work_path]:
                        if path and os.path.exists(path):
                            try:
                                shutil.rmtree(path)
                                print(f"成功删除目录: {path}")
                            except Exception as e:
                                print(f"删除目录失败: {path}, 错误: {str(e)}")
                
                sql_delete_tasks = "DELETE FROM tasks WHERE item_id = %s"
                cursor.execute(sql_delete_tasks, (item_id,))
                
                sql_delete_item = "DELETE FROM items WHERE item_id = %s"
                cursor.execute(sql_delete_item, (item_id,))
                
                conn.commit()
                return jsonify({"status": "success", "message": "项目删除成功"}), 200
        finally:
            conn.close()

    except Exception as e:
        print(f"删除项目失败: {str(e)}")
        return jsonify({"status": "error", "message": str(e)}), 500

@app.route('/api/delete_task', methods=['POST'])
def delete_task():
    try:
        item_id = request.headers.get('Item-Id')
        data = request.get_json()
        task_id = data.get('task_id', [])
        print(item_id,task_id)
        if not item_id or not task_id:
            return jsonify({"status": "error", "message": "参数不完整"}), 400

        conn = database.connect_db()
        try:
            with conn.cursor() as cursor:
                sql_check = "SELECT sh_id, result_path FROM tasks WHERE task_id = %s AND item_id = %s"
                cursor.execute(sql_check, (task_id, item_id))
                task = cursor.fetchone()
                
                if task:
                    sh_id, result_path = task
                    if sh_id and not check_slurm_job(str(sh_id)):
                        return jsonify({"status": "error", "message": "任务运行中，无法删除"}), 400
                    
                    paths_to_delete = []
                    if result_path:
                        try:
                            if isinstance(result_path, str):
                                try:
                                    path_dict = json.loads(result_path)
                                    if isinstance(path_dict, dict):
                                        paths_to_delete.extend(path_dict.values())
                                except json.JSONDecodeError:
                                    paths_to_delete.append(result_path)
                            elif isinstance(result_path, dict):
                                paths_to_delete.extend(result_path.values())
                        except Exception as e:
                            print(f"解析 result_path 失败: {str(e)}")
                            paths_to_delete = [result_path]

                    for path in paths_to_delete:
                        if isinstance(path, str) and os.path.exists(path):
                            try:
                                os.remove(path)
                                print(f"成功删除文件: {path}")
                            except Exception as e:
                                print(f"删除文件失败: {path}, 错误: {str(e)}")
                
                sql_delete = "DELETE FROM tasks WHERE task_id = %s AND item_id = %s"
                cursor.execute(sql_delete, (task_id, item_id))
                conn.commit()
                
                return jsonify({"status": "success", "message": "任务删除成功"}), 200
        finally:
            conn.close()

    except Exception as e:
        print(f"删除任务失败: {str(e)}")
        return jsonify({"status": "error", "message": str(e)}), 500
    

@app.route('/api/upload', methods=['POST'])
def upload_file():

    print("request.files:",request.files,'-'*20)
    if 'img' not in request.files:
        return jsonify(msg='找不到文件...'), 404

    file = request.files['img']
    if file.filename == '':
        return jsonify(msg='找不到文件...'), 404

    filename = secure_filename(file.filename)
    timestamp = int(time.time())
    new_filename = f"{timestamp}_{filename}"
    print(app.config['UPLOAD_FOLDER'])
    print(new_filename)
    final_filename = os.path.join(app.config['UPLOAD_FOLDER'], new_filename)
    file.save(final_filename)

    save_path = os.path.join(app.config['DATA_FOLDER'], new_filename[:-4])
    unzip(final_filename,save_path)

    save_path = os.path.join(save_path, filename[:-4])
    if check_files(save_path) == False:
        print(save_path)
        shutil.rmtree(save_path)
        return jsonify(msg='文件格式有误...'), 404
    
    figure_path = os.path.join(app.config['FIGURE_FOLDER'],new_filename[:-4])
    if not os.path.exists(figure_path):
        os.mkdir(figure_path)
        print(f"Folder {figure_path} created")
    else:
        print(f"Folder {figure_path} already exists")

    return jsonify(msg='ok', save_path= save_path, figure_path= figure_path), 200

@app.route('/api/upload_user_info', methods=['POST'])
def upload_user_info():
    try:
        item_id = request.headers.get('Item-Id')
        print("item_id", item_id)
        if not item_id:
            return jsonify({"msg": "未登录或登录已过期"}), 401

        work_path = database.get_item_workpath(item_id)[0]

        conn = database.connect_db()
        try:
            with conn.cursor() as cursor:
                sql = """SELECT save_path, species 
                        FROM items WHERE item_id = %s """
                cursor.execute(sql, (item_id,))
                items = cursor.fetchone()
                
                if items:
                    save_path = items[0]
                    species = items[1]
                else:
                    raise Exception("未找到对应的项目信息")
        finally:
            conn.close()
        figure_path = work_path

        num_nCount_RNA = request.form.get('num_nCount_RNA','')
        num_nFeature_RNA = request.form.get('num_nFeature_RNA','')
        num_percent_MT = request.form.get('num_percent_MT','25')  
        num_percent_HB = request.form.get('num_percent_HB','5')   

        timestamp = int(time.time() * 1000)
        Rpath = ""
        path_pdf = os.path.join(work_path, f"03.QC_{timestamp}.pdf")

        command = f'''Rscript {Rpath} -s "{save_path}" -f "{figure_path}" -p "{species}" -c "{num_nCount_RNA}" -r "{num_nFeature_RNA}" -m "{num_percent_MT}" -b "{num_percent_HB}" -o "{path_pdf}"'''
        print(command)

        param = {
            "item_id": item_id,
            "task_type": "parameter_setting",
            "parameter": {
                "num_nCount_RNA": num_nCount_RNA,
                "num_nFeature_RNA": num_nFeature_RNA,
                "num_percent_MT": num_percent_MT,
                "num_percent_HB": num_percent_HB,
            },
            "result_path": path_pdf,
        }

        sh_path = 'QCfinal.sh'
        job_id,_ = runsh_item(command, sh_path, param)
        print(f"任务{job_id}运行结束")

        return jsonify(msg='成功', task_id=job_id), 200

    except Exception as e:
        print(f"参数录入失败: {str(e)}")
        return jsonify({"msg": str(e)}), 500

@app.route('/api/get_figure_QC', methods=['GET'])
def get_figure_QC():
    try:
        item_id = request.headers.get('Item-Id')
        print("item_id", item_id)
        if not item_id:
            return jsonify({"msg": "未登录或登录已过期"}), 401

        work_path = database.get_item_workpath(item_id)[0]

        path_pdf1 = os.path.join(work_path, '01.QC.pdf')
        output_folder = os.path.join(work_path, 'figure')
        png_list1 = pdf_to_png(path_pdf1)

        png_list = png_list1
        print("结果很明显：", png_list)
        return jsonify({"png_paths": png_list}), 200
    except Exception as e:
        print(f"获取QC图失败: {str(e)}")
        return jsonify({"msg": str(e)}), 500
@app.route('/api/get_figure_QC_final', methods=['GET'])
def get_figure_QC_final():
    try:
        item_id = request.headers.get('Item-Id')
        print("item_id", item_id)
        if not item_id:
            return jsonify({"msg": "未登录或登录已过期"}), 401

        work_path = database.get_item_workpath(item_id)[0]

        task = database.get_latest_task_by_type(item_id=item_id, task_type="parameter_setting")
        if task:
            task_id, item_id, type, parameter, result_path, task_time, sh_id = task
            path_pdf = result_path
        else:
            path_pdf = os.path.join(work_path, '03.QC.pdf')
        png_list = pdf_to_png(path_pdf)
        print("最终结果很明显：", png_list)
        time.sleep(3)
        return jsonify({"png_paths": png_list}), 200
    except Exception as e:
        print(f"获取QC终图失败: {str(e)}")
        return jsonify({"msg": str(e)}), 500
    
@app.route('/api/get_SignificanceQC', methods=['POST'])
def get_SignificanceQC():
    try:
        data = request.get_json()
        selected_cells = data.get('selectedCells', [])
        
        if len(selected_cells) != 2:
            return jsonify({'error': '需要选择两个细胞类型'}), 400
        print(selected_cells)
            
        item_id = request.headers.get('Item-Id')
        print("item_id", item_id)
        if not item_id:
            return jsonify({"msg": "未登录或登录已过期"}), 401

        work_path = database.get_item_workpath(item_id)[0]
        cellid = selected_cells
        gene = data.get('selectedGene')

        path_pdf = os.path.join(pathf, 'QC_comparisons.pdf')
        image_path = pdf_to_svg(path_pdf)
        print("最终结果很明显：", image_path)
        time.sleep(2)

        return jsonify(msg='执行成功', success = True, imagePath = image_path), 200

    except Exception as e:
        return jsonify({'error': str(e)}), 500
    
@app.route('/api/get_figure_UMP', methods=['GET'])
def get_figure_UMP():
    try:
        item_id = request.headers.get('Item-Id')
        params = request.args.to_dict()
        task_id = params.get('params[taskId]', '')

        print("item_id",item_id)
        if not item_id:
            return jsonify({
                "status": "error",
                "message": "未登录或登录已过期"
            }), 401
        latest = False

        work_path = database.get_item_workpath(item_id)[0]
        if task_id:
            print("task_id",task_id)
            task = database.get_task_by_id(task_id)
            if task: 
                _, _, task_type, para, result_path, _, _ = task
                if task_type == 'Annotation': 
                    umap_path = result_path
                    ann_model = para["annotation_model"]
                else: latest = True
            else: latest = True
        else: latest = True
        if latest: 
            umap_path = os.path.join(work_path, 'umap_data.csv')
            ann_model = 'seurat_clusters'
            
        print("=="*50)
        print(umap_path)
        clustered_data,clustered_distribution,cluster_colors = read_csv(umap_path)
        print("duqvchenggong")
        clustered_data = clustered_data.to_json(orient='records') 
        return jsonify(msg='执行成功', clustered_data = clustered_data, clustered_distribution = clustered_distribution, cluster_colors=cluster_colors, ann_model = ann_model), 200
        
    except Exception as e:
        print(e)
        print(f"获取文件失败: {str(e)}")
        return jsonify({
            "status": "error",
            "message": str(e)
        }), 500


@app.route('/api/get_annotation_models', methods=['GET'])
def get_annotation_models():
    try:
        item_id = request.headers.get('Item-Id')
        print("item_id",item_id)
        if not item_id:
            return jsonify({
                "status": "error",
                "message": "未登录或登录已过期"
            }), 401

        work_path = database.get_item_workpath(item_id)[0]
        celltypist_models = list(models.models_description().model)
        children = []

        pkl_path = os.path.join(work_path,"reference_seurat_data","pkl")
        if os.path.exists(pkl_path):
            pkl_list = [f for f in os.listdir(pkl_path) if os.path.isfile(os.path.join(pkl_path, f))]
            celltypist_models = celltypist_models + pkl_list
        for model in celltypist_models:
            model = model[:-4]
            children.append({'value':model,
                            'label':model,})
        children.append({'value':"train",
                        'label':"训练你的模型",})
        celltypist = {
            'value': 'celltypist',
            'label': 'celltypist',
            'children':children}
        return jsonify(celltypist=celltypist), 200        
    except Exception as e:
        print(f"获取文件失败: {str(e)}")
        return jsonify({
            "status": "error",
            "message": str(e)
        }), 500


@app.route('/api/select_annotation_model', methods=['POST'])
def select_annotation_model():
    try:
        item_id = request.headers.get('Item-Id')
        print("item_id",item_id)
        if not item_id:
            return jsonify({
                "status": "error",
                "message": "未登录或登录已过期"
            }), 401

        work_path = database.get_item_workpath(item_id)[0]
        
        annotation_model = request.form.get('annotation_model')
        singleR_models = ['HumanPrimaryCellAtlasData','BlueprintEncodeData','DatabaseImmuneCellExpressionData','NovershternHematopoieticData','MonacoImmuneData','MouseRNAseqData','ImmGenData']
        timestamp = int(time.time() * 1000)
        print(annotation_model)
        task_id_new = ''
        if (annotation_model == "AICelltype"):
            print("选择AICelltype注释模型")
            categories = "AICelltype"
            def process_marker_data(file_path):
                """
                读取marker基因CSV文件，并根据cluster分组，将基因组织成列表的列表。

                Args:
                    file_path (str): marker基因CSV文件的路径。

                Returns:
                    list: 包含每个cluster的marker基因列表的列表。
                """
                try:
 
                    df = pd.read_csv(file_path)


                    if 'gene' not in df.columns or 'cluster' not in df.columns:
                        print("错误：CSV文件必须包含 'gene' 和 'cluster' 列。")
                        return []

                    gene_by_cluster = []
                    cluster_list = []
                    for cluster_name, group in df.groupby('cluster', sort=False):
                        genes = group['gene'].unique().tolist()
                        gene_by_cluster.append(genes)
                        cluster_list.append(cluster_name)
                                
                    return cluster_list, gene_by_cluster

                except FileNotFoundError:
                    print(f"错误：文件未找到，请检查路径是否正确：{file_path}")
                    return []
                except Exception as e:
                    print(f"处理文件时发生错误：{e}")
                    return []

            task = database.get_latest_task_by_type(item_id=item_id, task_type="Marker")
            if task:
                task_id, _, _, parameter, result_path, task_time, sh_id = task
                marker_gene_path = result_path["marker_data"]
            else:
                marker_gene_path = os.path.join(work_path, 'marker_data.csv')
            print("marker_gene_path:", marker_gene_path)

            result_cluster_list,result_gene_list = process_marker_data(marker_gene_path)

            gene_lt = result_gene_list
            tissue_name = database.get_item_notes(item_id)[0]

            cell_lt = aicelltype(tissue_name, gene_lt, model='qwen-max')
            print(cell_lt)
            anno = {
                "seurat_clusters":result_cluster_list,
                "annotation": cell_lt
            }
            print(anno)
            anno=DataFrame(anno)
            df6 = pd.merge(os.path.join(work_path,"umap_data.csv"),anno,how='right',on='seurat_clusters')
            df6.to_csv(f"{work_path}/umap_data_ai.csv",index=False)

        if (annotation_model in singleR_models):
            print("选择SingleR注释模型:",annotation_model)
            categories = "SingleR"
            
            umap_data_path = os.path.join(work_path,f'umap_data_{timestamp}.csv')

            param = {
                "item_id":item_id,
                "task_type":"Annotation",
                "parameter":{'annotation_model': annotation_model},
                "result_path":umap_data_path,
            }

            Rpath = "SingleRAnnotation.R"
            command = f'''Rscript {Rpath} -d "{work_path}" -m "{annotation_model}" -t "{timestamp}"'''
            print(command)
            sh_path = 'SingleR.sh'
            job_id,task_id_new = runsh_item(command, sh_path, param)
            print(f"任务{job_id}运行结束")
            print("完成注释")
            
        else:
            print("选择celltypist注释模型:",annotation_model)
            categories = "celltype"
            h5ad_path = os.path.join(work_path,'seurat_data.h5ad')
            
            if not os.path.isfile(h5ad_path):
                print("正在进行数据格式转换")
                Rpath = "convert_rds_to_h5ad.R"
                file = "seurat_data.rds"
                command = f'''Rscript {Rpath} -w "{work_path}" -f "{file}" -s "seurat_data.h5Seurat"'''
                print(command)
                sh_path = 'convert.sh'
                job_id = runsh(command, sh_path)
                print(f"任务{job_id}运行结束")
                print("转换完成")
                
            umap_data_path = copy_file_with_number(os.path.join(work_path,f'umap_data.csv'), timestamp)

            command = f'''celltypist_run.py --work_path "{work_path}" --annotation_model "{annotation_model}" --umap_data_path "{umap_data_path}"'''
            print(f"使用Python脚本进行注释: {command}")
            sh_path = 'celltypist.sh'

            job_id = runsh(command, sh_path)
            print(f"任务{job_id}运行结束")

            Rpath = "saveRDS.R"
            command = f'''Rscript {Rpath} -w "{work_path}" -f "{umap_data_path}"'''
            print(command)
            sh_path = 'saveRDS.sh'
            param = {
                "item_id":item_id,
                "task_type":"Annotation",
                "parameter":{'annotation_model': annotation_model},
                "result_path":umap_data_path,
            }

            job_id,task_id_new = runsh_item(command, sh_path, param)
            print(f"任务{job_id}运行结束")
            
        Rpath = "scale_diagram.R"
        seriesNames = "orig.ident"
        
        command = f'''Rscript {Rpath} -w "{work_path}" -s "{seriesNames}" -c "{categories}"'''
        print(command)
        sh_path = 'scale.sh'
        job_id = runsh(command, sh_path)
        print(f"任务{job_id}运行结束")
        print("完成注释")
        
        return jsonify({
                "task_id_new": task_id_new
            }), 200
    except Exception as e:
        print(f"注释失败: {str(e)}")
        return jsonify({
            "status": "error",
            "message": str(e)
        }), 500

@app.route('/api/large_model_annotation', methods=['POST'])
def large_model_annotation():
    try:
        start = datetime.now()

        item_id = request.headers.get('Item-Id')
        print("item_id",item_id)
        if not item_id:
            return jsonify({
                "status": "error",
                "message": "未登录或登录已过期"
            }), 401

        work_path = database.get_item_workpath(item_id)[0]
        data = request.get_json()
        task_id = data.get('taskId', '')

        print("task_id",task_id)
        task = database.get_task_by_id(task_id)
        if task: 
             _, _, task_type, param, result_path, _, _ = task
        else: 
            return jsonify({
                "status": "error",
                "message": "任务不存在"
            }), 401

        print(param)
        
        if "annotation_model" not in param.keys():
            return jsonify({
                "status": "error",
                "message": "未执行过注释任务，请先进行注释操作"
            }), 401
        annotation_model = param["annotation_model"]
        print(annotation_model)
        
        conn = database.connect_db()
        try:
            with conn.cursor() as cursor:
                sql = """
                    SELECT result_path 
                    FROM tasks 
                    WHERE item_id = %s AND type = %s 
                    AND JSON_EXTRACT(parameter, '$.annotation_model') = %s
                """
                cursor.execute(sql, (item_id, 'Marker', annotation_model))
                result = cursor.fetchone()
                if result:
                    result, = result
                    try:
                        if isinstance(result, str):
                            result = json.loads(result)
                    except json.JSONDecodeError:
                        pass  
                    print(result)
                    marker_gene_path = result["marker_data"]
                else:
                    return jsonify({
                        "status": "error",
                        "message": "未执行过计算marker基因的任务，请先计算marker基因"
                    }), 401
                
                cursor.execute(sql, (item_id, 'Annotation', annotation_model))
                result2 = cursor.fetchone()
                if result2:
                    result2, = result2
                    
                    print(result2)
                    umap_data_path = result
        finally:
            conn.close()

        print(marker_gene_path)
        def process_marker_data(file_path):
                """
                读取marker基因CSV文件，并根据cluster分组，将基因组织成列表的列表。

                Args:
                    file_path (str): marker基因CSV文件的路径。

                Returns:
                    list: 包含每个cluster的marker基因列表的列表。
                """
                try:
                    df = pd.read_csv(file_path)

                    if 'gene' not in df.columns or 'cluster' not in df.columns:
                        print("错误：CSV文件必须包含 'gene' 和 'cluster' 列。")
                        return []


                    gene_by_cluster = []
                    cluster_list = []
                    for cluster_name, group in df.groupby('cluster', sort=False):
                        genes = group['gene'].unique().tolist()
                        gene_by_cluster.append(genes)
                        cluster_list.append(cluster_name)
                                
                    return cluster_list, gene_by_cluster

                except FileNotFoundError:
                    print(f"错误：文件未找到，请检查路径是否正确：{file_path}")
                    return []
                except Exception as e:
                    print(f"处理文件时发生错误：{e}")
                    return []
        
        result_cluster_list,result_gene_list = process_marker_data(marker_gene_path)

        gene_lt = result_gene_list
        tissue_name = database.get_item_notes(item_id)[0]

        cell_lt = aicelltype(tissue_name, gene_lt, model='qwen-max')
        cell_lt = cell_lt = [x for x in cell_lt if x not in ('', '\n')]

        anno_new = []
        reason = []
        print("--"*50)
        for item in cell_lt:
            key, value = item.split(':')
            anno_new.append(key)
            reason.append(value)
        print(cell_lt)

        anno = {
            "origin":result_cluster_list,
            "large": anno_new,
            "reason": reason
        }
        anno=DataFrame(anno)
        anno = anno.to_dict(orient='records')
   
        end = datetime.now()
        elapsed_time = end - start
        print(f"代码运行时间: {elapsed_time}")

        return jsonify({"anno":anno}), 200
    except Exception as e:
        print(f"获取文件失败: {str(e)}")
        return jsonify({
            "status": "error",
            "message": str(e)
        }), 500

@app.route('/api/upload_reference', methods=['POST'])
def upload_reference():
    try:
        item_id = request.headers.get('Item-Id')
        print("item_id",item_id)
        if not item_id:
            return jsonify({
                "status": "error",
                "message": "未登录或登录已过期"
            }), 401

        work_path = database.get_item_workpath(item_id)[0]
        print("request.files:",request.files,'-'*20)
        if 'zip' not in request.files:
            return jsonify(msg='找不到文件...'), 404
        print(request.files)
        file = request.files['zip']
        if file.filename == '':
            return jsonify(msg='找不到文件...'), 404

        wkdir = os.path.join(work_path,"reference_seurat_data")
        upload_path = os.path.join(wkdir,'upload')
        file_path = os.path.join(wkdir,'data')
        pkl_path = os.path.join(wkdir,'pkl')

        if not os.path.exists(wkdir):
            print("这是用户在该任务首次上传训练数据，正在创建相关文件夹中……")
            os.makedirs(wkdir)
            print(f"已创建文件夹：{wkdir}")

            os.makedirs(upload_path)
            print(f"已创建文件夹：{upload_path}")

            os.makedirs(file_path)
            print(f"已创建文件夹：{file_path}")

            os.makedirs(pkl_path)
            print(f"已创建文件夹：{pkl_path}")

        filename = secure_filename(file.filename)
        upload_path = os.path.join(upload_path, filename) 
        print("开始保存...")
        file.save(upload_path)

        save_path = os.path.join(file_path, filename[:-4])
        unzip(upload_path,save_path)

        if os.path.exists(save_path):
            print("成功上传")
        
        return jsonify(msg='ok', save_path=save_path), 200
    except Exception as e:
        print(f"上传失败: {str(e)}")
        return jsonify({
            "status": "error",
            "message": str(e)
        }), 500

@app.route('/api/submit_reference_data', methods=['POST'])
def submit_reference_data():

    item_id = request.headers.get('Item-Id')
    print("item_id",item_id)
    if not item_id:
        return jsonify({
                "status": "error",
                "message": "未登录或登录已过期"
        }), 401

    work_path = database.get_item_workpath(item_id)[0]

    save_path = request.form.get('save_path')
    model_name = request.form.get('model_name')
    current_date = datetime.now().strftime("%Y-%m-%d")
    reference_dataset = os.path.join(save_path, "reference_seurat_data.h5ad")

    for root, dirs, files in os.walk(save_path):
        file = files[0]
        print("file:",file)
    Rpath = "convert_rds_to_h5ad.R"

    command = f'''Rscript {Rpath} -w "{save_path}" -f "{file}" -s "reference_seurat_data.h5Seurat"'''
    print(command)
    sh_path = 'convert.sh'
    job_id = runsh(command, sh_path)

    
    wait_count = 0
    while not os.path.exists(reference_dataset) and wait_count < 60:
        print(f"等待h5seurat文件生成: {reference_dataset}，已经等待{wait_count*5}秒")
        time.sleep(5)
        wait_count += 1
    print(f"任务{job_id}运行结束")

    adata = sc.read_h5ad(reference_dataset)
    custom_model = celltypist.train(X=adata, labels='Annotation',
        use_SGD=True, n_jobs=1,
        feature_selection=True,
        date=current_date)

    pkl_path = os.path.join(work_path,"reference_seurat_data","pkl")
    custom_model.write(f'{pkl_path}/{model_name}.pkl')
    return jsonify(msg=f'训练结束，模型名为{model_name}'), 200 


@app.route('/api/get_category_model', methods=['GET'])
def get_category_model():
    try:
        item_id = request.headers.get('Item-Id')
        print("item_id", item_id)
        if not item_id:
            return jsonify({
                "status": "error",
                "message": "未登录或登录已过期"
            }), 401

        work_path = database.get_item_workpath(item_id)[0]

        file_path = os.path.join(work_path,'category_models.txt')
        try:
            with open(file_path, 'r', encoding='utf-8') as file:
                content = file.read()
                print("文件的内容如下：")
                print(content)
        except FileNotFoundError:
            print(f"错误：文件 '{file_path}' 未找到。请检查文件路径是否正确。")
        except Exception as e:
            print(f"读取文件时发生错误：{e}")

        models = content.split('\n')[:-1]
        model_list = [{'value': model, 'label': model} 
                    for model in models]
        
        return jsonify({
            "model_list": model_list
        }), 200

    except Exception as e:
        print(f"获取基因列表失败: {str(e)}")
        return jsonify({
            "status": "error",
            "message": str(e)
        }), 500


@app.route('/api/get_scale_diagram_data', methods=['GET'])
def get_chart_data():
    try:
        item_id = request.headers.get('Item-Id')
        print("item_id", item_id)
        if not item_id:
            return jsonify({
                "status": "error",
                "message": "未登录或登录已过期"
            }), 401

        work_path = database.get_item_workpath(item_id)[0]
        

        params = request.args.to_dict()
        seriesNames = params.get('params[seriesNames]', 'seurat_clusters')
        categories = params.get('params[categoryMode]', 'orig.ident')


        models_file_path = os.path.join(work_path,'category_models.txt')
        try:
            with open(models_file_path, 'r', encoding='utf-8') as file:
                content = file.read()
                print("文件的内容如下：")
                print(content)
        except FileNotFoundError:
            print(f"错误：文件 '{file_path}' 未找到。请检查文件路径是否正确。")
        except Exception as e:
            print(f"读取文件时发生错误：{e}")
        if seriesNames=="seurat_clusters": 
            serie = seriesNames       
        else: serie = "celltype"

        print("categories:", categories)
        print("seriesNames:", serie)
        file_path = os.path.join(work_path, f"scale_diagram_{seriesNames}_{categories}.csv")
        if not os.path.exists(file_path):
            Rpath = "scale_diagram.R"
            command = f'''Rscript {Rpath} -w "{work_path}" -s "{serie}" -c "{categories}" -f "{file_path}"'''
            print(command)
            sh_path = 'scale.sh'
            job_id = runsh(command, sh_path)
            print(f"任务{job_id}运行结束")
        
            
        freq_df = pd.read_csv(file_path, index_col=0)

        rawData = freq_df.values.tolist()
        categories = freq_df.columns.tolist()
        seriesNames = freq_df.index.tolist()    
        print(freq_df)
        return jsonify({
            "status": "success",
            "rawData": rawData,
            "seriesNames": seriesNames,
            "categories": categories
        }), 200
        
    except Exception as e:
        print(f"获取文件失败: {str(e)}")
        return jsonify({
            "status": "error",
            "message": str(e)
        }), 500


@app.route('/api/calculate_marker_gene', methods=['POST'])
def calculate_marker_gene():
    try:
        item_id = request.headers.get('Item-Id')
        print("item_id",item_id)
        if not item_id:
            return jsonify({
                "status": "error",
                "message": "未登录或登录已过期"
            }), 401

        work_path = database.get_item_workpath(item_id)[0]
        annotation_model = "cluster"
        task = database.get_latest_task_by_type(item_id=item_id, task_type="Annotation")
        if task:
            task_id, item_id, type, parameter, result_path, task_time, sh_id = task
            annotation_model = parameter["annotation_model"]
            print("当前的注释模型：",annotation_model)
        timestamp = int(time.time() * 1000)
        all_marker_path = os.path.join(work_path,f'all_marker_{timestamp}.csv')
        marker_data_path = os.path.join(work_path,f'marker_data_{timestamp}.csv')
        all_genes_expression_path = os.path.join(work_path,f'all_genes_expression_{timestamp}.csv')

        Rpath = "FindMarkerGene.R"
        command = f'''Rscript {Rpath} -w "{work_path}" -t {timestamp}'''
        
        print(command)
        
        sh_path = 'marker_gene.sh'
        param = {
            "item_id": item_id,
            "task_type": "Marker",
            "parameter": {"annotation_model": annotation_model},
            "result_path": {"all_marker": all_marker_path, "marker_data": marker_data_path, "all_genes_expression": all_genes_expression_path},
        }

        job_id,_ = runsh_item(command, sh_path, param)
        time.sleep(5)
        print(f"任务{job_id}运行结束")
        
        return jsonify(msg='运行结束'), 200

    except Exception as e:
        print(f"计算marker基因失败: {str(e)}")
        return jsonify({
            "status": "error",
            "message": str(e)
        }), 500

@app.route('/api/get_figure_feature', methods=['POST'])
def get_figure_feature():
    try:
        item_id = request.headers.get('Item-Id')
        print("item_id",item_id)
        if not item_id:
            return jsonify({
                "status": "error",
                "message": "未登录或登录已过期"
            }), 401
        work_path = database.get_item_workpath(item_id)[0]

        gene = request.form.get('gene')
        group_model = 'celltype'
        timestamp = int(time.time() * 1000)
        pdf_path = os.path.join(work_path, f'marker_genes_{timestamp}.pdf')
        expression_path = os.path.join(work_path, f'gene_expression_{timestamp}.csv')

        task = database.get_latest_task_by_type(item_id=item_id, task_type="Annotation")
        if task:
            task_id, item_id, type, parameter, result_path, task_time, sh_id = task
            annotation_model = parameter["annotation_model"]
            print("当前的注释模型：",annotation_model)
        else:
            all_marker_path = os.path.join(work_path, 'all_marker.csv')

        Rpath = "FeaturePlot.R"
        command = f'''Rscript {Rpath} -w "{work_path}" -g "{gene}" -m "{group_model}" -t "{timestamp}"'''
        print(command)
        
        sh_path = 'feature_plot.sh'
        param = {
            "item_id": item_id,
            "task_type": "Feature",
            "parameter": {"group_model": group_model, "gene": gene, "annotation_model": annotation_model},
            "result_path": {"pdf": pdf_path, "expression": expression_path}
        }

        job_id,_ = runsh_item(command, sh_path, param)
        time.sleep(5)
        print(f"任务{job_id}运行结束")


        gene_expression_df = pd.read_csv(expression_path)
        gene_expression = list(np.around(gene_expression_df["gene_expression"], 3))

        img_list = pdf_to_png(pdf_path)
        
        return jsonify(msg='执行成功', gene_expression=gene_expression, img_list=img_list), 200

    except Exception as e:
        print(f"获取特征图失败: {str(e)}")
        return jsonify({
            "status": "error",
            "message": str(e)
        }), 500
@app.route('/api/fetch_feature_figure', methods=['POST'])
def fetch_feature_figure():
    try:
        item_id = request.headers.get('Item-Id')
        print("item_id",item_id)
        if not item_id:
            return jsonify({
                "status": "error",
                "message": "未登录或登录已过期"
            }), 401

        task_id = request.form.get('taskId', '')
        

        if task_id and task_id != 'undefined':
            task = database.get_task_by_id(task_id)
            if task:
                _, _, task_type, _, result_path, _, _ = task
                if task_type == "Feature":
                    expression_path = result_path["expression"]
                    pdf_path = result_path["pdf"]

        gene_expression_df = pd.read_csv(expression_path)
        gene_expression = list(np.around(gene_expression_df["gene_expression"], 3))
        print("=="*30)

        img_list = pdf_to_png(pdf_path)
        
        return jsonify(msg='执行成功', gene_expression=gene_expression, img_list=img_list), 200

    except Exception as e:
        print(f"获取特征图失败: {str(e)}")
        return jsonify({
            "status": "error",
            "message": str(e)
        }), 500


@app.route('/api/get_marker_genes', methods=['POST'])
def get_marker_genes():
    try:
        item_id = request.headers.get('Item-Id')
        if not item_id:
            return jsonify({
                "status": "error",
                "message": "未登录或登录已过期"
            }), 401
        work_path = database.get_item_workpath(item_id)[0]

        num = request.form.get('num', 5)
        task_id = request.form.get('taskId', '')
        latest = False
        file_path = ''
        if task_id:
            print("task_id", task_id)
            task = database.get_task_by_id(task_id)
            if task:
                _, _, task_type, _, result_path, _, _ = task
                if task_type == "Marker": file_path = result_path["all_marker"]
                else: latest = True
            else: latest = True
        else: latest = True
        if latest:
            if not item_id:
                return jsonify({
                    "status": "error",
                    "message": "未登录或登录已过期"
                }), 401
            print("item_id", item_id)

            task = database.get_latest_task_by_type(item_id=item_id, task_type="Marker")
            if task:
                task_id, item_id, type, parameter, result_path, task_time, sh_id = task
                file_path = result_path["all_marker"]
            else:
                file_path = os.path.join(work_path, 'all_marker.csv')
                
        print("file_path:", file_path,"num:", num)
        data_markers = pd.read_csv(file_path)
        
        top_markers_py = (data_markers.groupby('cluster')
                         .apply(lambda x: x.nlargest(int(num), 'avg_log2FC'))
                         .reset_index(drop=True))
        
        gene_list = top_markers_py.to_dict(orient='records')

        return jsonify({"gene_list":gene_list}), 200

    except Exception as e:
        print(f"获取marker基因失败: {str(e)}")
        return jsonify({
            "status": "error",
            "message": str(e)
        }), 500


@app.route('/api/get_figure_dotplot', methods=['POST'])
def get_figure_dotplot():
    try:
        item_id = request.headers.get('Item-Id')
        if not item_id:
            return jsonify({
                "status": "error",
                "message": "未登录或登录已过期"
            }), 401

        work_path = database.get_item_workpath(item_id)[0]

        gene_num = int(request.form.get('num', 5))
        task_id = request.form.get('taskId', '')
        print("task_id的值是：",task_id, "gene_num的值是：", gene_num)
        latest = False

        if task_id:
            print("taskId",task_id)
            task = database.get_task_by_id(task_id)
            if task:
                _, _, task_type, _, result_path, _, _ = task
                if task_type == "Marker":
                    all_marker_path = result_path["all_marker"]
                    all_genes_expression_path = result_path["all_genes_expression"]
                else: latest = True
            else: latest = True
        else: latest = True
        if latest:
            print("item_id",item_id)
            task = database.get_latest_task_by_type(item_id=item_id, task_type="Marker")
            if task:
                task_id, item_id, type, parameter, result_path, task_time, sh_id = task
                all_marker_path = result_path["all_marker"]
                all_genes_expression_path = result_path["all_genes_expression"]
            else:
                all_marker_path = os.path.join(work_path, 'all_marker.csv')
                all_genes_expression_path = os.path.join(work_path, 'all_genes_expression.csv')

        marker_data = pd.read_csv(all_marker_path)
        top_genes = marker_data.groupby('cluster').apply(
            lambda x: x.nlargest(gene_num, 'avg_log2FC')
        ).reset_index(drop=True)['gene'].unique()


        expression_data = pd.read_csv(all_genes_expression_path)
        dotplot_data = expression_data[expression_data['Gene'].isin(top_genes)]
        
        result = dotplot_data.to_json(orient='records')
        
        return jsonify(msg='执行成功', dotplot_data=result), 200

    except Exception as e:
        print(f"获取气泡图失败: {str(e)}")
        return jsonify({
            "status": "error", 
            "message": str(e)
        }), 500


@app.route('/api/get_figure_barplot', methods=['POST'])
def get_figure_barplot():
    try:
        item_id = request.headers.get('Item-Id')
        print("item_id",item_id)
        if not item_id:
            return jsonify({
                "status": "error",
                "message": "未登录或登录已过期"
            }), 401

        work_path = database.get_item_workpath(item_id)[0]

        choose_cluster = request.form.get('choose_cluster')
        task_id = request.form.get('taskId', '')
        print("task_id的值是：",task_id, "cluster的值是：", choose_cluster)
        timestamp = int(time.time() * 1000)

        annotation_model = "cluster"
        task = database.get_latest_task_by_type(item_id=item_id, task_type="Marker")
        if task:
            task_id, item_id, type, parameter, result_path, task_time, sh_id = task
            all_marker_path = result_path["all_marker"]
            annotation_model = parameter["annotation_model"]
            print("当前的注释模型：",annotation_model)
        else:
            all_marker_path = os.path.join(work_path, 'all_marker.csv')
        go_data_path = os.path.join(work_path, f'go_data_{timestamp}.csv')

        Rpath = "Go.R"
        command = f'''Rscript {Rpath} -w "{work_path}" -c "{choose_cluster}" -m "{all_marker_path}" -t "{timestamp}"'''
        print(command)
                
        sh_path = 'go.sh'
        param = {
                "item_id": item_id,
                "task_type": "Go",
                "parameter": {"cluster": choose_cluster, "annotation_model":annotation_model},
                "result_path": go_data_path
        }

        job_id,_ = runsh_item(command, sh_path, param)
        print(f"任务{job_id}运行结束")
        time.sleep(3)

        bar_data = pd.read_csv(go_data_path)
        bar_data = bar_data.sort_values(by='Count', ascending=False).head(10)
        bar_data = bar_data.values.tolist()
        
        return jsonify(msg='执行成功', bar_data=bar_data), 200

    except Exception as e:
        print(f"获取GO图失败: {str(e)}")
        return jsonify({
            "status": "error",
            "message": str(e)
        }), 500

@app.route('/api/fetch_go_figure', methods=['POST'])
def fetch_go_figure():
    try:
        item_id = request.headers.get('Item-Id')
        print("item_id",item_id)
        if not item_id:
            return jsonify({
                "status": "error",
                "message": "未登录或登录已过期"
            }), 401

        task_id = request.form.get('taskId', '')
        latest = False
        if task_id:
            print("taskId",task_id)

            task = database.get_task_by_id(task_id)
            if task:
                _, _, task_type, _, result_path, _, _ = task
                if task_type == "Go": go_data_path = result_path
                else: latest = True
            else: latest = True
        else: latest = True
        if latest:
            task = database.get_latest_task_by_type(item_id=item_id, task_type="Go")
            if task:
                task_id, item_id, type, parameter, result_path, task_time, sh_id = task
                go_data_path = result_path
            else:
                return jsonify({
                        "status": "error",
                        "message": "暂无历史任务"
                    }), 404


        bar_data = pd.read_csv(go_data_path)
        bar_data = bar_data.sort_values(by='Count', ascending=False).head(10)
        bar_data = bar_data.values.tolist()
        
        return jsonify(msg='执行成功', bar_data=bar_data), 200

    except Exception as e:
        print(f"获取GO图失败: {str(e)}")
        return jsonify({
            "status": "error",
            "message": str(e)
        }), 500

@app.route('/api/perform_cellchat_analysis', methods=['POST'])
def perform_cellchat_analysis():
    try:
        item_id = request.headers.get('Item-Id')
        print("item_id",item_id)
        if not item_id:
            return jsonify({
                "status": "error",
                "message": "未登录或登录已过期"
            }), 401
        
        annotation_model = "cluster"
        work_path = database.get_item_workpath(item_id)[0]
        task = database.get_latest_task_by_type(item_id=item_id, task_type="Annotation")
        if task:
            task_id, item_id, type, parameter, result_path, task_time, sh_id = task
            umap_data_path = result_path
            annotation_model = parameter["annotation_model"]
            print("当前的注释模型：",annotation_model)
        else:
            umap_data_path = os.path.join(work_path, 'umap_data.csv')


        timestamp = int(time.time() * 1000)
        
        Rpath = "Cellchat.R"
        command = f'''Rscript {Rpath} -w "{work_path}" -a "{umap_data_path}" -t "{timestamp}"'''
        print(command)
            
        sh_path = 'cellchat.sh'
        result_paths = {
            "circle_data_weight": os.path.join(work_path, f"weight_{timestamp}.json"),
            "circle_data_count": os.path.join(work_path, f"count_{timestamp}.json"),
            "heatmap_weight": os.path.join(work_path, f"cellchat_heatmap_weight_{timestamp}.pdf"),
            "heatmap_count": os.path.join(work_path, f"cellchat_heatmap_count_{timestamp}.pdf")
        }
        
        param = {
            "item_id": item_id,
            "task_type": "CellChat",
            "parameter": {annotation_model:"annotation_model"},
            "result_path": result_paths
        }

        job_id,_ = runsh_item(command, sh_path, param)
        print(f"任务{job_id}运行结束")
        time.sleep(3)
            
        return jsonify({
            'success': True,
            'message': '分析完成'
        })
        
    except Exception as e:
        print(f"细胞通讯分析失败: {str(e)}")
        return jsonify({
            'success': False,
            'message': str(e)
        }), 500
    
@app.route('/api/get_circle_chart_data')
def get_circle_chart_data():
    try:
        item_id = request.headers.get('Item-Id')
        if not item_id:
            return jsonify({
                "status": "error",
                "message": "未登录或登录已过期"
            }), 401

        work_path = database.get_item_workpath(item_id)[0]
        params = request.args.to_dict()
        weight_type = params.get('params[weightType]', '')
        task_id = params.get('params[taskId]', '')
        latest = False
        file_path = ''
        if task_id:
            print("task_id",task_id)
            task = database.get_task_by_id(task_id)
            if task: 
                _, _, task_type, _, result_path, _, _ = task
                if task_type == 'CellChat': file_path = result_path[f"circle_data_{weight_type}"] 
                else: latest = True
            else: latest = True
        else: latest = True
        if latest: 
            task = database.get_latest_task_by_type(item_id=item_id, task_type="CellChat")
            print("圆图数据信息：",task)
            if task: 
                _, _, _, _, result_path, _, _ = task
                file_path = result_path[f"circle_data_{weight_type}"]
            else: 
                return jsonify({
                        "status": "error",
                        "message": "暂无历史任务"
                    }), 404

        with open(file_path, 'r') as f:
            data = json.load(f)
        return jsonify(data)
    except Exception as e:
        print(f"获取圆形图数据失败: {str(e)}")
        return jsonify({'error': str(e)}), 500

@app.route('/api/get_cellchat_heatmaps', methods=['GET'])
def get_cellchat_heatmaps():
    try:
        item_id = request.headers.get('Item-Id')
        if not item_id:
            return jsonify({
                "status": "error",
                "message": "未登录或登录已过期"
            }), 401

        work_path = database.get_item_workpath(item_id)[0]
        params = request.args.to_dict()
        weight_type = params.get('params[weightType]', '')
        task_id = params.get('params[taskId]', '')
        latest = False
        file_path = ''
        if task_id:
            print("task_id",task_id)

            task = database.get_task_by_id(task_id)
            if task: 
                _, _, task_type, _, result_path, _, _ = task
                if task_type == 'CellChat': file_path = result_path[f"heatmap_{weight_type}"] 
                else:
                    latest = True
            else: 
                latest = True
        else: latest = True
        if latest: 
            task = database.get_latest_task_by_type(item_id=item_id, task_type="CellChat")
            if task: 
                _, _, _, _, result_path, _, _ = task
                print(task)
                file_path = result_path[f"heatmap_{weight_type}"]
            else: 
                return jsonify({
                        "status": "error",
                        "message": "暂无历史任务"
                    }), 404
    
        img_list = pdf_to_png(file_path)
        print("最终结果很明显：", img_list)
        time.sleep(2)
        
        return jsonify(msg='执行成功', img_list=img_list), 200

    except Exception as e:
        print(f"获取热图失败: {str(e)}")
        return jsonify({
            "status": "error",
            "message": str(e)
        }), 500


@app.route('/api/submit_selections', methods=['POST'])
def submit_selections():
    try:
        item_id = request.headers.get('Item-Id')
        print("item_id", item_id)
        if not item_id:
            return jsonify({
                "status": "error",
                "message": "未登录或登录已过期"
            }), 401

        work_path = database.get_item_workpath(item_id)[0]
        data = request.get_json()
        first_selection = data.get('firstSelection', [])
        second_selection = data.get('secondSelection', [])

        print("First Selection:", first_selection)
        print("Second Selection:", second_selection)
        Rpath = "Infercnv.R"
        input_path = "seurat_data.rds"
        SubsetList = ','.join(first_selection)
        ReferenceList = ','.join(second_selection)
        timestamp = int(time.time() * 1000)
        command = f'''Rscript {Rpath} -w "{work_path}" -i "{input_path}" -s "{SubsetList}" -r "{ReferenceList}"'''
        sh_path = 'script.sh'

        result_path = os.path.join(work_path, 'infercnv.png')
        destination_directory = "src/assets/result"
        dst_path = os.path.join(destination_directory, f'infercnv_{timestamp}.png')

        annotation_model = "cluster"
        task = database.get_latest_task_by_type(item_id=item_id, task_type="Annotation")
        if task:
            task_id, item_id, type, parameter, _, task_time, sh_id = task
            annotation_model = parameter["annotation_model"]
            print("当前的注释模型：",annotation_model)
        param = {
            "item_id": item_id,
            "task_type": "InferCNV",
            "parameter": {
                "first_selection": first_selection,
                "second_selection": second_selection,
                "annotation_model": annotation_model
            },
            "result_path": f'infercnv_{timestamp}.png',
        }
        

        job_id,_ = runsh_item(command, sh_path, param)
        time.sleep(5)
        print(f"任务{job_id}运行结束")
        print(os.path.exists(result_path))
        print(result_path)
        if not os.path.exists(result_path):
            return jsonify({
                "status": "error",
                "message": "结果图片未生成"
            }), 500
        
        print(os.path.exists(result_path))
        os.makedirs(destination_directory, exist_ok=True)
        shutil.copy2(result_path, dst_path)
        os.remove(result_path)
        print("原图片已删除")

        return jsonify({
            'message': 'Selections received successfully',
            'imgPath': f'infercnv_{timestamp}.png',
        }), 200

    except Exception as e:
        print(f"提交infercnv选择失败: {str(e)}")
        return jsonify({
            "status": "error",
            "message": str(e)
        }), 500
    
@app.route('/api/get_infercnv_figure', methods=['GET'])
def get_infercnv_figure():
    try:
        item_id = request.headers.get('Item-Id')
        if not item_id:
            return jsonify({
                "status": "error",
                "message": "未登录或登录已过期"
            }), 401

        work_path = database.get_item_workpath(item_id)[0]
        params = request.args.to_dict()
        task_id = params.get('params[taskId]', '')
        latest = False

        if task_id:
            print("task_id",task_id)
            task = database.get_task_by_id(task_id)
            print(task)
            if task: 
                _, _, task_type, _, result_path, _, _ = task
                if task_type == 'InferCNV': file_path = result_path 
                else:
                    latest = True
            else: 
                latest = True
        else: latest = True 
        if latest: 
            task = database.get_latest_task_by_type(item_id=item_id, task_type="InferCNV")
            if task: 
                _, _, _, _, result_path, _, _ = task
                file_path = result_path
        print("infercnv图片路径：", file_path)
        return jsonify({
            'message': 'Selections received successfully',
            'imgPath': file_path,
        }), 200

    except Exception as e:
        print(f"获取InferCNV图失败: {str(e)}")
        return jsonify({
            "status": "error",
            "message": str(e)
        }), 500


@app.route('/api/submit_pseudotime_cells', methods=['POST'])
def submit_pseudotime_cells():
    try:
        item_id = request.headers.get('Item-Id')
        print("item_id", item_id)
        if not item_id:
            return jsonify({
                "status": "error",
                "message": "未登录或登录已过期"
            }), 401

        work_path = database.get_item_workpath(item_id)[0]


        annotation_model = "cluster"

        task = database.get_latest_task_by_type(item_id=item_id, task_type="Marker")
        if task:
            task_id, item_id, type, parameter, result_path, task_time, sh_id = task
            all_marker_path = result_path["all_marker"]
            annotation_model = parameter["annotation_model"]
            print("当前的注释模型：",annotation_model)
        else:
            all_marker_path = os.path.join(work_path, 'all_marker.csv')
        data = request.get_json()
        selection = data.get('selection', [])
        print("Selection:", selection)

        timestamp = int(time.time() * 1000)
        Rpath = " monocle.R"
        SelectionList = ','.join(selection)
        command = f'''Rscript {Rpath} -w "{work_path}" -s "{SelectionList}" -m "{all_marker_path}" -t "{timestamp}"'''
        print(command)

        path_pdf_tra = os.path.join(work_path, f"001.trajectory_{timestamp}.pdf")
        path_pdf_anno = os.path.join(work_path, f"001.annot_trajectory_{timestamp}.pdf")
        param = {
            "item_id": item_id,
            "task_type": "Pseudotime",
            "parameter": {"SelectionList": SelectionList, "annotation_model": annotation_model},
            "result_path": {
                "path_pdf_tra": path_pdf_tra,
                "path_pdf_anno": path_pdf_anno,
            },
        }

        sh_path = 'monocle.sh'
        job_id,_ = runsh_item(command, sh_path, param)
        time.sleep(5)
        print(f"任务已完成，ID: {job_id}")

        image_path_tra = pdf_to_png(path_pdf_tra)
        image_path_anno = pdf_to_png(path_pdf_anno)
        print("最终结果很明显：", image_path_tra, image_path_anno)
        time.sleep(2)

        return jsonify({
            "status": "success",
            "imgPaths_tra": image_path_tra,
            "imgPaths_anno": image_path_anno,
        }), 200

    except Exception as e:
        print(f"发生错误: {str(e)}")
        return jsonify({
            "status": "error",
            "message": str(e)
        }), 500


@app.route('/api/get_pseudo_genes', methods=['GET'])
def get_pseudo_genes():
    try:
        item_id = request.headers.get('Item-Id')
        print("item_id", item_id)
        if not item_id:
            return jsonify({
                "status": "error",
                "message": "未登录或登录已过期"
            }), 401

        work_path = database.get_item_workpath(item_id)[0]
        task = database.get_latest_task_by_type(item_id=item_id, task_type="Pseudotime")
        if task:
            task_id, item_id, type, parameter, result_path, task_time, sh_id = task
            path_pdf_anno = result_path["path_pdf_anno"]
            timestamp = path_pdf_anno.split('_')[-1][:-4]
        else:
            return jsonify({
                "status": "error",
                "message": "请先执行拟时序分析"
            }), 401
        
        file_path = os.path.join(work_path,f'HSMM_gene_{timestamp}.txt')
        try:
            with open(file_path, 'r', encoding='utf-8') as file:
                content = file.read()
                print("文件的内容如下：")
        except FileNotFoundError:
            print(f"错误：文件 '{file_path}' 未找到。请检查文件路径是否正确。")
        except Exception as e:
            print(f"读取文件时发生错误：{e}")

        all_genes = content.split("\n")[:-1]

        gene_list = [{'value': i, 'label': i} for i in all_genes]

        return jsonify({"gene_list": gene_list}), 200

    except Exception as e:
        print(f"获取pseudo基因列表失败: {str(e)}")
        return jsonify({
            "status": "error",
            "message": str(e)
        }), 500

@app.route('/api/selection_pseudotime_genes', methods=['POST'])
def selection_pseudotime_genes():
    try:
        item_id = request.headers.get('Item-Id')
        print("item_id", item_id)
        if not item_id:
            return jsonify({
                "status": "error",
                "message": "未登录或登录已过期"
            }), 401

        work_path = database.get_item_workpath(item_id)[0]
        data = request.get_json()
        selected_items = data.get('selected_items', [])
        print("selected_items:", selected_items)

        timestamp = int(time.time() * 1000)
        Rpath = " monocleGene.R"
        file_path = os.path.join(work_path, "HSMM_object.rds")
        pdf_path = os.path.join(work_path, f"003.pseudotime_plot_{timestamp}.pdf")
        SelectionList = ','.join(selected_items)

        command = f'''Rscript {Rpath} -i "{file_path}" -p "{pdf_path}" -s "{SelectionList}"'''
        sh_path = 'monoclegene.sh'
        annotation_model = "cluster"
        task = database.get_latest_task_by_type(item_id=item_id, task_type="Annotation")
        if task:
            task_id, item_id, type, parameter, result_path, task_time, sh_id = task
            annotation_model = parameter["annotation_model"]
        param = {
            "item_id": item_id,
            "task_type": "Pseudotime_gene",
            "parameter": {
                "SelectionList": SelectionList,
                "annotation_model": annotation_model
            },
            "result_path": pdf_path,
        }

        job_id,_ = runsh_item(command, sh_path, param)
        time.sleep(5)


        image_path_gene = pdf_to_png(pdf_path)
        time.sleep(2)
        
        return jsonify({
            "status": "success",
            "imgPaths_gene": image_path_gene
        }), 200

    except Exception as e:
        print(f"发生错误: {str(e)}")
        return jsonify({
            "status": "error",
            "message": str(e)
        }), 500
@app.route('/api/query_latest_result_pseudotime', methods=['GET'])
def query_latest_result_pseudotime():
    try:
        item_id = request.headers.get('Item-Id')
        if not item_id:
            return jsonify({
                "status": "error",
                "message": "未登录或登录已过期"
            }), 401
        params = request.args.to_dict()
        task_id = params.get('params[taskId]', '')
        path_pdf_tra = None
        path_pdf_anno = None
        path_pdf_gene = None
        latest = False

        if task_id:

            task = database.get_task_by_id(task_id)
            if task: 
                _, _, task_type, _, result_path, _, _ = task
                if task_type == 'Pseudotime':
                    path_pdf_tra = result_path.get("path_pdf_tra")
                    path_pdf_anno = result_path.get("path_pdf_anno")                
                elif task_type == 'Pseudotime_gene':
                    path_pdf_gene = result_path
                else:
                    latest = True
            else: 
                latest = True
        else: latest = True
        if latest: 
            task = database.get_latest_task_by_type(item_id=item_id, task_type="Pseudotime")
            if task:
                task_id, item_id, type, parameter, result_path, task_time, sh_id = task
                path_pdf_tra = result_path.get("path_pdf_tra")
                path_pdf_anno = result_path.get("path_pdf_anno")
            else:
                return jsonify({
                        "status": "error",
                        "message": "暂无历史任务"
                    }), 404
            task_gene = database.get_latest_task_by_type(item_id=item_id, task_type="Pseudotime_gene")
            if task_gene:
                _, _, _, _, result_path_gene, _, _ = task_gene
                path_pdf_gene = result_path_gene


        image_path_tra = pdf_to_png(path_pdf_tra)
        image_path_anno = pdf_to_png(path_pdf_anno)
        image_path_gene = pdf_to_png(path_pdf_gene)
        time.sleep(2)
        
        return jsonify({
            "status": "success",
            "imgPaths_tra": image_path_tra,
            "imgPaths_anno": image_path_anno,
            "imgPaths_gene": image_path_gene
        }), 200

    except Exception as e:
        print(f"发生错误: {str(e)}")
        return jsonify({
            "status": "error",
            "message": str(e)
        }), 500

@app.route('/api/selection_tcga_gene', methods=['POST'])
def selection_tcga_gene():
    try:
        item_id = request.headers.get('Item-Id')
        print("item_id", item_id)
        if not item_id:
            return jsonify({
                "status": "error",
                "message": "未登录或登录已过期"
            }), 401

        data = request.get_json()
        query = data.get('query', '')
        work_path = database.get_item_workpath(item_id)[0]
        latest = False

        path_pdf_boxplot = None
        path_pdf_survival = None
        path_pdf_meth = None
        path_pdf_TIMER = None
        table_csv = None
        if query:
            print("正在查询TCGA的结果")
            task_id = data.get('taskId', '')
            if task_id:
                print("TCGA的task_id：",task_id)

                task = database.get_task_by_id(task_id)
                if task: 
                    _, _, task_type, _, result_path, _, _ = task
                    if task_type == 'TCGA':
                        path_pdf_boxplot = result_path.get("path_pdf_boxplot")
                        path_pdf_survival = result_path.get("path_pdf_survival")
                        path_pdf_meth = result_path.get("path_pdf_meth")
                        path_pdf_TIMER = result_path.get("path_pdf_TIMER")
                        table_csv = result_path.get("table_csv")               
                    else: latest = True
                else: 
                    latest = True
            else: latest = True
            if latest:
                task = database.get_latest_task_by_type(item_id=item_id, task_type="TCGA")
                if task:
                    _, _, _, parameter, result_path, _, _ = task
                    path_pdf_boxplot = result_path.get("path_pdf_boxplot")
                    path_pdf_survival = result_path.get("path_pdf_survival")
                    path_pdf_meth = result_path.get("path_pdf_meth")
                    path_pdf_TIMER = result_path.get("path_pdf_TIMER")
                    table_csv = result_path.get("table_csv")
                else:
                    return jsonify({
                        "status": "error",
                    }), 404
        else:
            print("正在训练TCGA结果")
            gene = data.get('selected_items', '')
            timestamp = int(time.time() * 1000)
            Rpath = [
                " TCGA/code/4.TCGA_mRNA_survival.R",
                " TCGA/code/5.TCGA_mRNA_meth.R",
                " TCGA/code/6.TCGA_mRNA_TIMER2.0.R",
                " TCGA/code/7.TCGA_bayes.R"
            ]

            wd_path = work_path

            path_pdf_boxplot = os.path.join(wd_path, 'TCGA', f'2.TCGA_mRNA_boxplot_{timestamp}.pdf')
            path_pdf_survival = os.path.join(wd_path, 'TCGA', f'4.TCGA_mRNA_survival_{timestamp}.pdf')
            path_pdf_meth = os.path.join(wd_path, 'TCGA', f'5.TCGA_mRNA_meth_{timestamp}.pdf')
            path_pdf_TIMER = os.path.join(wd_path, 'TCGA', f'6.TCGA_mRNA_TIMER2.0_{timestamp}.pdf')
            table_csv = os.path.join(wd_path, f'TCGA/res/{gene}_bayes_imm&TPM_{timestamp}.csv')

            rpath = " TCGA/code/2.TCGA_mRNA_boxplot.R"
            command = f'''Rscript {rpath} -w " TCGA" -s "{wd_path}"  -g "{gene}" -t "{timestamp}"'''
            sh_path = 'tcga.sh'

            annotation_model = "cluster"
            task = database.get_latest_task_by_type(item_id=item_id, task_type="Annotation")
            if task:
                task_id, item_id, type, parameter, result_path, task_time, sh_id = task
                annotation_model = parameter["annotation_model"]
                print("当前的注释模型：",annotation_model)
            param = {
                    "item_id": item_id,
                    "task_type": "TCGA",
                    "parameter": {
                        "gene": gene,
                        "annotation_model": annotation_model
                    },
                    "result_path": {
                        "path_pdf_boxplot": path_pdf_boxplot,
                        "path_pdf_survival": path_pdf_survival,
                        "path_pdf_meth": path_pdf_meth,
                        "path_pdf_TIMER": path_pdf_TIMER,
                        "table_csv": table_csv
                    },
            }
            job_id,_ = runsh_item(command, sh_path, param)
            time.sleep(5)
            for idx, rpath in enumerate(Rpath):
                command = f'''Rscript {rpath} -w " TCGA" -s "{wd_path}"  -g "{gene}" -t "{timestamp}"'''
                sh_path = 'tcga.sh'
                job_id = runsh_item_update(command, sh_path, job_id)
                time.sleep(5)

        img_boxplot = pdf_to_png(path_pdf_boxplot)
        img_survival = pdf_to_png(path_pdf_survival)
        img_meth = pdf_to_png(path_pdf_meth)
        img_TIMER = pdf_to_png(path_pdf_TIMER)
 
        if os.path.exists(table_csv):
            tabledata = pd.read_csv(table_csv)
            tabledata = tabledata.fillna("")
            tabledata = tabledata.to_dict(orient='records')
        else:
            tabledata = []

        time.sleep(2.5)
        return jsonify({
            "status": "success",
            "img_boxplot": img_boxplot,
            "img_survival": img_survival,
            "img_meth": img_meth,
            "img_TIMER": img_TIMER,
            "tabledata": tabledata
        }), 200

    except Exception as e:

        msg = str(e)
        if "no such file" in msg or "No such file" in msg:
            import re
            match = re.search(r"no such file: '([^']+)'", msg, re.IGNORECASE)
            if match:
                missing_file = match.group(1)
                missing_file = missing_file.split("/")[-1]
                return jsonify({
                    "status": "error",
                    "message": f"缺少文件{missing_file}，请尝试重新提交"
                }), 500
            else:
                return jsonify({
                    "status": "error",
                    "message": f"缺少文件，具体信息：{msg}，请尝试重新提交"
                }), 500
        return jsonify({
            "status": "error",
            "message": msg
        }), 500


@app.route('/api/selection_tcga_genes', methods=['POST'])
def selection_tcga_genes():
    try:
        item_id = request.headers.get('Item-Id')
        if not item_id:
            return jsonify({
                "status": "error",
                "message": "未登录或登录已过期"
            }), 401

        data = request.get_json()
        query = data.get('query', '')
        work_path = database.get_item_workpath(item_id)[0]
        latest = False
        path_pdf_heatmap = None
        forest_dir = None

        if query:
            task_id = data.get('taskId', '')
            if task_id:
                task = database.get_task_by_id(task_id)
                if task: 
                    _, _, task_type, _, result_path, _, _ = task
                    if task_type == 'TCGA_genes':
                        path_pdf_heatmap = result_path.get("path_pdf_heatmap")
                        forest_dir = result_path.get("forest_dir")              
                    else: latest = True
                else: 
                    latest = True
            else: latest = True
            if latest:
                task = database.get_latest_task_by_type(item_id=item_id, task_type="TCGA_genes")
                if task:
                    _, _, _, parameter, result_path, _, _ = task
                    path_pdf_heatmap = result_path.get("path_pdf_heatmap")
                    forest_dir = result_path.get("forest_dir")
                else:
                    return jsonify({
                        "status": "error",
                        "message": "暂无历史任务"
                    }), 404
        else:
            gene = data.get('selected_items', [])
            gene_str = ','.join(gene)
            timestamp = int(time.time() * 1000)
            Rpath = " TCGA/code/8.TCGA_genelist_heatmap.R"
            wd_path = work_path
            forest_dir = os.path.join(wd_path, f'TCGA/plot/forest_list_{timestamp}')
            path_pdf_heatmap = os.path.join(wd_path, 'TCGA', f'08.heatmap_{timestamp}.pdf')

            if os.path.exists(forest_dir):
                shutil.rmtree(forest_dir)
            os.makedirs(forest_dir, exist_ok=True)

            command = f'''Rscript {Rpath} -w " TCGA" -s "{wd_path}" -g "{gene_str}" -t "{timestamp}"'''

            sh_path = 'tcga.sh'
            annotation_model = "cluster"
            task = database.get_latest_task_by_type(item_id=item_id, task_type="Annotation")
            if task:
                task_id, item_id, type, parameter, result_path, task_time, sh_id = task
                annotation_model = parameter["annotation_model"]
            param = {
                "item_id": item_id,
                "task_type": "TCGA_genes",
                "parameter": {
                    "gene": gene,
                    "annotation_model": annotation_model
                },
                "result_path": {
                    "path_pdf_heatmap": path_pdf_heatmap,
                    "forest_dir": forest_dir
                },
            }
            job_id,_ = runsh_item(command, sh_path, param)

            time.sleep(5)

        img_heatmap = pdf_to_png(path_pdf_heatmap)
        img_forest = []
        if os.path.exists(forest_dir):
            for item in os.listdir(forest_dir):
                item_path = os.path.join(forest_dir, item)
                img = pdf_to_png(item_path)
                img_forest += img


        time.sleep(2.5)
        return jsonify({
            "status": "success",
            "img_heatmap": img_heatmap,
            "img_forest": img_forest,
        }), 200

    except Exception as e:
        return jsonify({
            "status": "error",
            "message": str(e)
        }), 500
 
@app.route('/api/tcga_resize', methods=['POST'])
def tcga_resize():
    try:
        task_id = request.form.get('task_id')
        outlier_size = request.form.get('outlier_size')         
        frame_size = request.form.get('frame_size')
        box_width = request.form.get('box_width')
        axis_textsize = request.form.get('axis_textsize')
        base_size = request.form.get('base_size')
        ylimits = request.form.get('ylimits')
        print("task_id:",task_id)
        print(outlier_size,frame_size,box_width,axis_textsize,base_size,ylimits)

        item_id = request.headers.get('Item-Id')
        print("item_id", item_id)
        if not item_id:
            return jsonify({
                "status": "error",
                "message": "未登录或登录已过期"
            }), 401

        work_path = database.get_item_workpath(item_id)[0]

        task = database.get_task_by_id(task_id)
        _, _, _, para, result_path, _, _ = task


        Rpath = " TCGA/code/resize.R"
        file_path = result_path["path_pdf_boxplot"]
        gene = para["gene"]

        command = f'''Rscript {Rpath} -w "{work_path}" -g "{gene}" -i "{file_path}" -o "{outlier_size}" -s "{frame_size}" -b "{box_width}" -a "{axis_textsize}" -z "{base_size}" -y "{ylimits}"'''
        sh_path = 'run.sh'
        job_id = runsh(command, sh_path)
        time.sleep(5)
        print(f"任务{job_id}运行结束")
        return jsonify(msg='成功'), 200


    except Exception as e:

        msg = str(e)
        if "no such file" in msg or "No such file" in msg:

            import re
            match = re.search(r"no such file: '([^']+)'", msg, re.IGNORECASE)
            if match:
                missing_file = match.group(1)
                missing_file = missing_file.split("/")[-1]
                return jsonify({
                    "status": "error",
                    "message": f"缺少文件{missing_file}，请尝试重新提交"
                }), 500
            else:
                return jsonify({
                    "status": "error",
                    "message": f"缺少文件，具体信息：{msg}，请尝试重新提交"
                }), 500

        return jsonify({
            "status": "error",
            "message": msg
        }), 500



@app.route('/api/get_all_tasks', methods=['GET'])
def get_all_tasks():
    conn = sqlite3.connect('tasks.db')
    c = conn.cursor()
    c.execute("SELECT task_id FROM tasks")
    tasks = c.fetchall()
    conn.close()

    task_ids = [task[0] for task in tasks]
    return jsonify(task_ids=task_ids), 200

@app.route('/api/get_gene_list', methods=['GET'])
def get_gene_list():
    try:
        item_id = request.headers.get('Item-Id')
        print("item_id", item_id)
        if not item_id:
            return jsonify({
                "status": "error",
                "message": "未登录或登录已过期"
            }), 401

        work_path = database.get_item_workpath(item_id)[0]

        task = database.get_latest_task_by_type(item_id=item_id, task_type="Marker")
        if task:
            task_id, item_id, type, parameter, result_path, task_time, sh_id = task
            all_marker_path = result_path["all_marker"]
        else:
            all_marker_path = os.path.join(work_path, 'all_marker.csv')
        
        data_markers = pd.read_csv(all_marker_path)
        

        top_markers = (data_markers.groupby('cluster')
                      .apply(lambda x: x.nlargest(20, 'avg_log2FC'))
                      .reset_index(drop=True))

        gene_list = [{'value': gene, 'label': gene} 
                    for gene in top_markers['gene'].unique()]
        

        gene_list_all = [{'value': gene, 'label': gene} 
                        for gene in data_markers['gene'].unique()]
        
        return jsonify({
            "gene_list": gene_list,
            "gene_list_all": gene_list_all
        }), 200

    except Exception as e:
        print(f"获取基因列表失败: {str(e)}")
        return jsonify({
            "status": "error",
            "message": str(e)
        }), 500

@app.route('/api/get_cell_list', methods=['GET'])
def get_cell_list():
    try:
        item_id = request.headers.get('Item-Id')
        print("item_id", item_id)
        if not item_id:
            return jsonify({
                "status": "error",
                "message": "未登录或登录已过期"
            }), 401

        work_path = database.get_item_workpath(item_id)[0]

        task = database.get_latest_task_by_type(item_id=item_id, task_type="Annotation")
        if task:
            task_id, item_id, type, parameter, result_path, task_time, sh_id = task
            umap_data_path = result_path
            print(parameter)
        else:
            umap_data_path = os.path.join(work_path, 'umap_data.csv')
        celllist = pd.read_csv(umap_data_path)
        celllist = celllist['annotation'].unique().tolist()
        return jsonify({'celllist': celllist}), 200
    except Exception as e:
        print(e)
        return jsonify({'error': str(e)}), 500

@app.route('/api/login', methods=['POST'])
def login():
    try:
        user_name = request.form.get('user_name')
        password = request.form.get('password')
        

        conn = pymysql.connect(host=host, port=3306, user='sunpc', password='spc20020629', database='scseq')
        cursor = conn.cursor()
        

        cursor.execute('SELECT id, name FROM users WHERE name = %s AND password = %s', (user_name, password))
        result = cursor.fetchone()
        
        if result:
            print(f"'user_id': {result[0]},'user_name':{ result[1]}")
            return jsonify({
                'status': 'success',
                'message': '登录成功',
                'user_id': result[0],
                'user_name': result[1]
            }), 200
        else:
            return jsonify({
                'status': 'error',
                'message': '用户名或密码错误'
            }), 401
            
    except Exception as e:
        return jsonify({
            'status': 'error',
            'message': str(e)
        }), 500
        
    finally:
        if 'cursor' in locals():
            cursor.close()
        if 'conn' in locals():
            conn.close()

@app.route('/api/register', methods=['POST'])
def register():
    try:
        user_name = request.form.get('user_name')
        password = request.form.get('password')
        phone = request.form.get('phone')
        
        conn = pymysql.connect(host=host, port=3306, user='sunpc', password='spc20020629', database='scseq')
        cursor = conn.cursor()
        
        cursor.execute('SELECT id FROM users WHERE name = %s', (user_name,))
        if cursor.fetchone():
            return jsonify({
                'status': 'error',
                'message': '用户名已存在'
            }), 400
            
        cursor.execute('SELECT id FROM users WHERE phone = %s', (phone,))
        if cursor.fetchone():
            return jsonify({
                'status': 'error',
                'message': '该手机号已被注册'
            }), 400
        
        cursor.execute(
            'INSERT INTO users (name, password, phone) VALUES (%s, %s, %s)',
            (user_name, password, phone)
        )
        conn.commit()
        print("注册成功")
        
        return jsonify({
            'status': 'success',
            'message': '注册成功'
        }), 200
            
    except Exception as e:
        return jsonify({
            'status': 'error',
            'message': str(e)
        }), 500
        
    finally:
        if 'cursor' in locals():
            cursor.close()
        if 'conn' in locals():
            conn.close()



@app.route('/api/check_user_info', methods=['GET'])
def check_user_info():
    try:
        user_id = request.headers.get('User-Id')
        if not user_id:
            return jsonify({
                "success": False,
                "message": "未检测到用户信息"
            }), 401
            
        return jsonify({
            "success": True,
            "user_id": user_id,
            "message": "用户信息获取成功"
        })
        
    except Exception as e:
        return jsonify({
            "success": False,
            "message": str(e)
        }), 500
    
if __name__ == '__main__':
    ipv4s = socket.gethostbyname_ex(socket.gethostname())[2]
    print(ipv4s[0])
    app.run(debug=True,host=ipv4s[0],port=5001)
