<template>
  <el-card>
    <el-row>
      <el-col :span="12">
        <div class="grid-content">
          <span class="demonstration">QC图</span>
          <div class="demo-image__lazy">
            <div class="insert"></div>
            <el-scrollbar ref="el_scrollbar" style="height: 500px; overflow-y: auto">
              <el-image
              v-for="path in imgPaths"
                :key="path"
                :src="require('../assets/result/' + path)"
                @load="onImageLoad"
                @error="onImageError"
                style="width: 400px; height: 400px; margin: 10px"
                fit="contain"
                lazy
                :scroll-container="scrollContainer"
                :preview-src-list="[require('../assets/result/' + path)]"
              >
                <template #toolbar="{ actions, reset }">
                  <el-icon @click="actions('zoomOut')">
                    <ZoomOut />
                  </el-icon>
                  <el-icon
                    @click="actions('zoomIn', { enableTransition: false, zoomRate: 2 })"
                  >
                    <ZoomIn />
                  </el-icon>
                  <el-icon
                    @click="
                      actions('clockwise', { rotateDeg: 180, enableTransition: false })
                    "
                  >
                    <RefreshRight />
                  </el-icon>
                  <el-icon @click="actions('anticlockwise')">
                    <RefreshLeft />
                  </el-icon>
                  <el-icon @click="reset">
                    <Refresh />
                  </el-icon>
                  <el-icon @click="downloadImage(require('../assets/result/' + path))">
                    <Download />
                  </el-icon>
                </template>
                <template #error>
                  <div class="image-slot">
                    <el-icon><icon-picture /></el-icon>
                  </div>
                </template>
              </el-image>
            </el-scrollbar>
            <div class="insert"></div>
          </div>
        </div>
      </el-col>

      <!-- <el-divider direction="vertical" class="height: 100%;"></el-divider> -->

      <el-col :span="12">
        <div class="grid-content">
          <span class="demonstration"
            >参数设置<el-button
              type=""
              text
              size="small"
              @click="jumpToTutorial('parameter-settings')"
            >
              <el-icon><QuestionFilled /></el-icon> </el-button
          ></span>
          <el-form :model="numberValidateForm" ref="numberValidateForm">
            <el-form-item
              label="nCount_RNA>"
              prop="num_nCount_RNA"
              :rules="[
                { required: true, message: 'nCount_RNA不能为空' },
                { type: 'number', message: 'nCount_RNA必须为数字值' },
              ]"
            >
              <el-input
                v-model.number="numberValidateForm.num_nCount_RNA"
                autocomplete="off"
              ></el-input>
            </el-form-item>
            <el-form-item
              label="nFeature_RNA>"
              prop="num_nFeature_RNA"
              :rules="[
                { required: true, message: 'nFeature_RNA不能为空' },
                { type: 'number', message: 'nFeature_RNA必须为数字值' },
              ]"
            >
              <el-input
                v-model.number="numberValidateForm.num_nFeature_RNA"
                autocomplete="off"
              ></el-input>
            </el-form-item>

            <el-form-item>
              <el-form-item
                label="percent_MT<"
                prop="num_percent_MT"
                :rules="numberValidateForm.if_num_percent_MT ? [
                  { required: true, message: 'percent_MT不能为空' },
                  { type: 'number', message: 'percent_MT必须为数字值' },
                ] : []"
              >
                <el-input
                  v-model.number="numberValidateForm.num_percent_MT"
                  :disabled="!numberValidateForm.if_num_percent_MT"
                  autocomplete="off"
                ></el-input>
              </el-form-item>
              <div class="insert_in_form"></div>

              <el-form-item label="是否启用">
                <el-switch v-model="numberValidateForm.if_num_percent_MT"></el-switch>
              </el-form-item>
            </el-form-item>
            <el-form-item>
              <el-form-item
                label="percent_HB<"
                prop="num_percent_HB"
                :rules="numberValidateForm.if_num_percent_HB ? [
                  { required: true, message: 'percent_HB不能为空' },
                  { type: 'number', message: 'percent_HB必须为数字值' },
                ] : []"
              >
                <el-input
                  v-model.number="numberValidateForm.num_percent_HB"
                  :disabled="!numberValidateForm.if_num_percent_HB"
                  autocomplete="off"
                ></el-input>
              </el-form-item>
              <div class="insert_in_form"></div>

              <el-form-item label="是否启用">
                <el-switch v-model="numberValidateForm.if_num_percent_HB"></el-switch>
              </el-form-item>
            </el-form-item>
            <!-- <el-form-item label="ncol">
              <el-input v-model="numberValidateForm.ncol"></el-input>
            </el-form-item>

            <el-form-item label="pt_size">
              <el-input v-model="numberValidateForm.pt_size"></el-input>
            </el-form-item>

            <el-form-item label="scale_factor">
              <el-input v-model="numberValidateForm.scale_factor"></el-input>
            </el-form-item>

            <el-form-item label="nfeatures">
              <el-input v-model="numberValidateForm.nfeatures"></el-input>
            </el-form-item>

            <el-form-item label="dims_use">
              <el-input v-model="numberValidateForm.dims_use"></el-input>
            </el-form-item>

            <el-form-item label="lambda">
              <el-input v-model="numberValidateForm.lambda"></el-input>
            </el-form-item>

            <el-form-item label="resolution">
              <el-input v-model="numberValidateForm.resolution"></el-input>
            </el-form-item> -->

            <el-form-item label="聚类方法">
              <el-select
                v-model="numberValidateForm.clustering_method"
                placeholder="聚类方法"
              >
                <el-option label="方法一" value="1"></el-option>
                <el-option label="方法二" value="2"></el-option>
              </el-select>
            </el-form-item>

            <!-- <el-form-item label="聚类参数">
              <el-input v-model="numberValidateForm.resolution"></el-input>
            </el-form-item> -->

            <el-form-item>
              <el-row :gutter="20">
                <el-col :span="6"></el-col>

                <el-col :span="6"
                  ><el-button type="primary" @click="submitForm('numberValidateForm')"
                    >提交</el-button
                  ></el-col
                >
                <el-col :span="6"></el-col>
                <el-col :span="6"
                  ><el-button @click="resetForm('numberValidateForm')"
                    >重置</el-button
                  ></el-col
                >
              </el-row>
            </el-form-item>
          </el-form>
        </div>
      </el-col>
    </el-row>
  </el-card>
</template>

<script>
import axios from "axios";
import { ElMessage, ElMessageBox } from "element-plus";
import { Eleme } from "@element-plus/icons-vue";
import { Picture as IconPicture } from "@element-plus/icons-vue";
import { get, post } from "@/request/api";
import { UploadFilled } from "@element-plus/icons-vue";

import {
  Download,
  Refresh,
  RefreshLeft,
  RefreshRight,
  ZoomIn,
  ZoomOut,
} from "@element-plus/icons-vue";

export default {
  components: {
    Eleme,
    IconPicture,
    UploadFilled,
    Download,
    Refresh,
    RefreshLeft,
    RefreshRight,
    ZoomIn,
    ZoomOut,
  },
  data() {
    return {
      imgPaths: [],

      numberValidateForm: {
        num_nCount_RNA: "",
        num_nFeature_RNA: "",
        num_percent_MT: "",
        num_percent_HB: "",
        ncol: "",
        pt_size: "",
        scale_factor: "",
        nfeatures: "",
        dims_use: "",
        lambda: "",
        resolution: "",
        if_num_percent_MT: false,
        if_num_percent_HB: false,
        clustering_method: "",
      },
      task_id: this.$myContent.task_id,
      // scrollContainer: HTMLElement
    };
  },
  methods: {
    submitForm(formName) {
      const formData = new FormData();
      formData.append("num_nCount_RNA", this.numberValidateForm.num_nCount_RNA);
      formData.append("num_nFeature_RNA", this.numberValidateForm.num_nFeature_RNA);
      if (this.numberValidateForm.if_num_percent_HB) {
        formData.append("num_percent_HB", this.numberValidateForm.num_percent_HB);
      }
      if (this.numberValidateForm.if_num_percent_MT) {
        formData.append("num_percent_MT", this.numberValidateForm.num_percent_MT);
      }

      this.$refs[formName].validate((valid) => {
        if (valid) {
          alert("submit!");
          post("upload_user_info", formData)
            .then((response) => {
              ElMessage({
                message: "参数信息已录入",
                type: "success",
              });
              console.log(response.data);
              // 可根据返回的 task_id 做跳转或状态管理
              this.$router.push("/jsjg");
            });
        } else {
          console.log("error submit!!");
          return false;
        }
      });
    },
    resetForm(formName) {
      this.$refs[formName].resetFields();
    },

    onImageError() {
      console.log("Failed to load image");
      this.errorMessage = "Failed to load image";
      get(`get_figure_QC`)
        .then((res) => {
          console.log("res.data", res.data);
          console.log("res.data.png_paths", res.data.png_paths);
          this.imgPaths = res.data.png_paths;
        })
        .catch((error) => {
          console.log("error", error);
          if (error.response) {
            console.log("请求失败:", error.response.status);
            console.log("错误响应数据:", error.response.data);
            ElMessage.error(error.response.data.msg);

            // console.error("qwer", error);
          } else {
            console.error("请求错误:", error.message);
          }
        });
      // 你可以在这里处理加载错误，比如显示错误消息或尝试加载备用图片
    },

    downloadImage(url) {
      if (!url) return;
      console.log(url);

      // 从URL中提取文件名
      const filename = url.split("/").pop() || "image_" + Date.now() + ".jpeg";

      // 使用原始URL直接下载
      fetch(url)
        .then((response) => {
          if (!response.ok) {
            throw new Error("Network response was not ok");
          }
          return response.blob();
        })
        .then((blob) => {
          const blobUrl = window.URL.createObjectURL(blob);
          const link = document.createElement("a");
          link.href = blobUrl;
          link.download = filename;
          document.body.appendChild(link);
          link.click();
          document.body.removeChild(link);
          window.URL.revokeObjectURL(blobUrl);
        })
        .catch((error) => {
          console.error("Download failed:", error);
          ElMessage.error("下载失败：" + error.message);
        });
    },

    autoExecuteFunction() {
      if (this.imgPaths.length === 0) {
        this.onImageError();
      }
    },

    // 添加跳转方法
    jumpToTutorial(section) {
      // 修改路径为根路径
      window.open(`/tutorial.html#${section}`, "_blank");
    },
  },
  mounted() {
    // this.scrollContainer = this.$refs.el_scrollbar.wrap;  // 将 el-scrollbar 的 wrap 对象找出来，指定给 scroll-container
    // 当组件被挂载到 DOM 上之后，调用 autoExecuteFunction
    this.autoExecuteFunction();
  },
};
</script>

<style scoped>
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

.demonstration {
  display: block;
  color: var(--el-text-color-secondary);
  font-size: 14px;
  margin-bottom: 20px;
  text-align: center;
}

.image-slot {
  display: flex;
  justify-content: center;
  align-items: center;
  width: 100%;
  height: 100%;
  background: var(--el-fill-color-light);
  color: var(--el-text-color-secondary);
  font-size: 30px;
}

.el-row {
  margin-bottom: 20px;

  &:last-child {
    margin-bottom: 0;
  }
}

.el-col {
  border-radius: 4px;
}

.bg-purple-dark {
  background: #99a9bf;
}

.bg-purple {
  background: #fbfbfb;
}

.bg-purple-light {
  background: #e5e9f2;
}

.grid-content {
  border-radius: 4px;
  min-height: 36px;
  margin: 0 auto;
  /* 居中显示 */
}

.row-bg {
  padding: 10px 0;
  background-color: #f9fafc;
}

.demo-image__lazy {
  display: block;
  /* min-height: 200px; */
  margin-bottom: 200px;
  height: 100%;
  justify-content: space-between;
  margin-top: 0px;
  display: flex;
}

.insert {
  width: 20%;
}

.insert_in_form {
  width: 10%;
}
</style>
