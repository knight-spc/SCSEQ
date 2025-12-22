import pandas as pd
import matplotlib.pyplot as plt
import numpy as np



def rgb_to_hex(rgb):
    r, g, b = rgb
    
    r_int = int(r * 255)
    g_int = int(g * 255)
    b_int = int(b * 255)
    r = max(0, min(255, r_int))
    g = max(0, min(255, g_int))
    b = max(0, min(255, b_int))
    hex_string = '#{:02X}{:02X}{:02X}'.format(r, g, b)
    return hex_string

def generate_colors(num_colors, cmap_name='Spectral_r'):
    """
    生成给定个数的颜色。

    参数:
    num_colors (int): 要生成的颜色个数。
    cmap_name (str): matplotlib的颜色映射名称（默认为'viridis'）。gist_rainbow,rainbow_r,Set3,Paired,Spectral_r,jet_r好用。Reds,cool,warm可以用来做气泡图

    返回:
    list: 包含RGB元组的颜色列表。
    """
    # 获取颜色映射
    cmap = plt.cm.get_cmap(cmap_name)
    
    # 生成颜色列表（RGBA格式，但这里只取RGB）
    colors = [cmap(i / (num_colors - 1))[:3] for i in range(num_colors)]
    
    # 如果需要将颜色转换为十六进制字符串，可以这样做（但通常不需要）
    colors_hex = [rgb_to_hex(color) for color in colors]
    
    # 返回颜色列表（RGB格式）
    return colors, colors_hex



def read_csv(file_path):
    data = pd.read_csv(file_path)
    # 数据和标签
    x = data['umap_data.umap_1']
    y = data['umap_data.umap_2']
    labels = data['seurat_clusters']
    data_distribution = {
        'xmin': x.min(),
        'xmax': x.max(),
        'ymin': y.min(),
        'ymax': y.max(),
    }

    # 使用函数生成对应个数的颜色
    colors_rgb, colors_hex = generate_colors(labels.nunique())
    print(len(colors_hex))
    print(colors_hex)
    # color_list = [colors_hex[i-1] for i in labels]
    return data,data_distribution,colors_hex


