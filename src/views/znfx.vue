<template>
  <el-card>
    <div class="form-container">
      <el-form :model="form" label-width="90px" style="width: 60%;">

        <div class="demo-image__error">
          <div class="block">
            <span class="demonstration">UMAP图</span>
            <el-image :src="imageSrc_res" fit="contain" @load="onImageLoad" @error="onImageError">
              <template #error>
                <div class="image-slot">
                  <el-icon><icon-picture /></el-icon>
                </div>
              </template>
            </el-image>
          </div>
        </div>

      </el-form>

    </div>
    

  </el-card>
</template>
  
  
<script>
import axios from 'axios';
import { ElMessage, ElMessageBox } from 'element-plus'
import { Eleme } from '@element-plus/icons-vue'
import { Picture as IconPicture } from '@element-plus/icons-vue'

import { h } from 'vue'

export default {
  components: {
    Eleme,
    IconPicture
  },

  data() {
    return {
      task_id: this.$myContent.task_id,
      imageSrc_res: '',
      if_load: false,
    }
  },
  methods: {
    diagnosis() {
      //创建formdata对象
      // let that = this
      // // 确认是否缺少信息
      // if (this.task_id==' ') {
      //   console.error('请选择一个用户');
      //   ElMessage.error('请选择一个用户')
      //   return;
      // }
      this.if_load = true


      // 开始诊断
      // axios.get(`/api/eval?task_id=${this.task_id}`)
      axios.get(`/api/eval?task_id=${this.task_id}`)
        .then((res) => {
          console.log("查看图像中", res)
          // 我想要分割图
          axios.get(`/api/segment?task_id=${this.task_id}`, { responseType: 'arraybuffer' }).
            then((res) => {
              console.log("返回值：", res)
              let binaryData = res.data
              let blob = new Blob([binaryData], { type: 'image/jpeg' })
              var url = window.URL.createObjectURL(blob)
              // that.imageSrc = url
              this.imageSrc_res = url
            })

        })
        .catch(error => {
          console.log('error', error)
          if (error.response) {
            console.log('请求失败:', error.response.status);
            console.log('错误响应数据:', error.response.data);
            ElMessage.error(error.response.data.msg)

            // console.error("qwer", error);
          } else {
            console.error('请求错误:', error.message);
          }
        })

    },


    fetchUserTasks() {
      axios
        .get('/api/get_all_tasks')
        .then((response) => {
          this.userTasks = response.data.task_ids;
        })
        .catch((error) => {
          console.error('获取用户列表失败', error);
        })
    },
    fetchUserImg() {
      // 我想要原始图片
      axios.get(`/api/origin?task_id=${this.task_id}`, { responseType: 'arraybuffer' }).
        then((res) => {
          console.log("yuanshi", res)
          let binaryData = res.data
          let blob = new Blob([binaryData], { type: 'image/jpeg' })
          var url = window.URL.createObjectURL(blob)
          // that.imageSrc = url
          this.imageSrc_res = url
        })
    },
    onImageLoad() {
      console.log('Image loaded successfully');
      // 你可以在这里执行其他操作，比如隐藏加载提示  
      
    },
    onImageError() {
      console.log('Failed to load image');
      this.errorMessage = 'Failed to load image';
      axios.get(`/api/get_figure_UMP`, { responseType: 'arraybuffer' }).
        then((res) => {
          console.log("yuanshi", res.data)
          let binaryData = res.data
          let blob = new Blob([binaryData])
          var url = window.URL.createObjectURL(blob)
          // that.imageSrc = url
          this.imageSrc_res = url
        })
        .catch(error => {
          console.log('error', error)
          if (error.response) {
            console.log('请求失败:', error.response.status);
            console.log('错误响应数据:', error.response.data);
            ElMessage.error(error.response.data.msg)

            // console.error("qwer", error);
          } else {
            console.error('请求错误:', error.message);
          }
        })
      // 你可以在这里处理加载错误，比如显示错误消息或尝试加载备用图片  
    },


  },
}
</script>

<style scoped>
.form-container {
  display: flex;
  justify-content: center;
  align-items: center;
  font-size: 24px;
}

.el-button .custom-loading .circular {
  margin-right: 6px;
  width: 18px;
  height: 18px;
  animation: loading-rotate 2s linear infinite;
}

.el-button .custom-loading .circular .path {
  animation: loading-dash 1.5s ease-in-out infinite;
  stroke-dasharray: 90, 150;
  stroke-dashoffset: 0;
  stroke-width: 2;
  stroke: var(--el-button-text-color);
  stroke-linecap: round;
}

.demo-image__error .block {
  padding: 30px 0;
  text-align: center;
  /* border-right: solid 1px var(--el-border-color); */
  display: inline-block;
  width: 49%;
  box-sizing: border-box;
  vertical-align: top;
}

.demo-image__error .blockright {
  padding: 30px 0;
  text-align: center;
  display: inline-block;
  width: 49%;
  box-sizing: border-box;
  vertical-align: top;
}

.demo-image__error .demonstration {
  display: block;
  color: var(--el-text-color-secondary);
  font-size: 14px;
  margin-bottom: 20px;
}

.demo-image__error .el-image {
  padding: 0 5px;
  max-width: 300px;
  max-height: 200px;
  width: 100%;
  height: 200px;
}

.demo-image__error .image-slot {
  display: flex;
  justify-content: center;
  align-items: center;
  width: 100%;
  height: 100%;
  background: var(--el-fill-color-light);
  color: var(--el-text-color-secondary);
  font-size: 30px;
}

.demo-image__error .image-slot .el-icon {
  font-size: 30px;
}</style>
  
  