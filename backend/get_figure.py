import fitz  # PyMuPDF 库
import os
import re
from PIL import Image
import shutil
import os
import svgwrite

pdf_path = r'' 
output_folder = r'' 

svg_list = [] 
png_list = []

def pdf_to_svg(pdf_path):
    output_folder = "src/assets/result"

    task = pdf_path.split("/")[-2]+'_'+pdf_path.split('/')[-1].split('.')[0]
       
    pdf_document = fitz.open(pdf_path)

    svg_list = []

    for page_num in range(pdf_document.page_count):
        page = pdf_document.load_page(page_num)

        svg_data = page.get_svg_image()

        svg_filename = task + "_" + f"page_{page_num + 1}.svg"
        svg_list.append(svg_filename)
        svg_save_path = os.path.join("src/assets/result", svg_filename)
        with open(svg_save_path, 'w', encoding='utf-8') as svg_file:
            svg_file.write(svg_data)

    pdf_document.close()
    print(f"PDF 已成功转换为 SVG, 并保存到文件夹: {output_folder}")
    print(svg_list)
    return svg_list

def pdf_to_png(pdf_path):
    output_folder = "src/assets/result"
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
    if (pdf_path == None): return

    floder_path = pdf_path.split("/")[-2]
    if(floder_path=="TCGA"):
        floder_path = pdf_path.split("/")[-3]+'_'+floder_path
    file_path = pdf_path.split('/')[-1]
    task = floder_path +'_'+file_path.split('.')[0]+'_'+file_path.split('.')[1]
      
    pdf_document = fitz.open(pdf_path)
    
    png_list = []
    
    for page_num in range(pdf_document.page_count):
        page = pdf_document.load_page(page_num)
        
        pix = page.get_pixmap(dpi=300, alpha=False)  
        
        img = Image.frombytes("RGB", [pix.width, pix.height], pix.samples)
        
        background = Image.new("RGB", img.size, (255, 255, 255))
        background.paste(img, mask=None)
        
        png_filename = task + "_" + f"page_{page_num + 1}.png"
        
        png_list.append(png_filename)
        
        background.save(os.path.join("src/assets/result", png_filename))
    
    pdf_document.close()
    print(f"PDF 已成功转换为 PNG, 并保存到文件夹: {output_folder}")
    print(png_list)
    return png_list

