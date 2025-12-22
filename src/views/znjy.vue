<template>
  <div class="photo-gallery">
    <!-- 在顶部添加测试按钮 -->
    <el-button type="primary" @click="checkUserInfo" style="margin-bottom: 20px;">
      检查用户信息
    </el-button>
    <div class="gallery-row">
      <div class="photo-item">
    <h3>图片列表<el-button type="" text size="small" @click="jumpToTutorial('interface-overview')">
      <el-icon><QuestionFilled /></el-icon>
    </el-button></h3>
    
        <el-scrollbar
          ref="el_scrollbar"
          style="height: 300px; overflow-y: auto; width: 100%"
        >
          <el-image
            v-for="(url, index) in imageList"
            :key="index"
            :src="url"
            style="height: 450px; margin: 10px auto; display: block"
            fit="contain"
            lazy
            :scroll-container="scrollContainer"
            :preview-src-list="[url]"
          >
            <template #toolbar="{ actions, reset }">
              <el-icon @click="actions('zoomOut')"><ZoomOut /></el-icon>
              <el-icon
                @click="actions('zoomIn', { enableTransition: false, zoomRate: 2 })"
                ><ZoomIn
              /></el-icon>
              <el-icon
                @click="actions('clockwise', { rotateDeg: 180, enableTransition: false })"
                ><RefreshRight
              /></el-icon>
              <el-icon @click="actions('anticlockwise')"><RefreshLeft /></el-icon>
              <el-icon @click="reset"><Refresh /></el-icon>
              <el-icon @click="downloadImage(url)"><Download /></el-icon>
            </template>
            <template #error>
              <div class="image-slot">
                <el-icon><icon-picture /></el-icon>
              </div>
            </template>
          </el-image>
        </el-scrollbar>
      </div>
      <div class="photo-item">
        <h3>图片 2</h3>
        <el-image :src="url2" :preview-src-list="[url2]" fit="contain">
          <!-- 工具栏模板相同 -->
        </el-image>
      </div>
    </div>
    <div class="gallery-row">
      <div class="photo-item">
        <h3>图片 3</h3>
        <el-image :src="url3" :preview-src-list="[url3]" fit="contain">
          <!-- 工具栏模板相同 -->
        </el-image>
      </div>
      <div class="photo-item">
        <h3>图片 4</h3>
        <el-image :src="url4" :preview-src-list="[url4]" fit="contain">
          <!-- 工具栏模板相同 -->
        </el-image>
      </div>
    </div>
  </div>
</template>

<script>
import { ElMessage } from "element-plus";
import { Delete, Edit, Search, Share, Upload, QuestionFilled } from '@element-plus/icons-vue'
import { get } from '@/request/api'  // 改用 get 方法

import {
  Download,
  Refresh,
  RefreshLeft,
  RefreshRight,
  ZoomIn,
  ZoomOut,
  Picture as IconPicture,
} from "@element-plus/icons-vue";

export default {
  components: {
    Download,
    Refresh,
    RefreshLeft,
    RefreshRight,
    ZoomIn,
    ZoomOut,
    IconPicture,
  },
  data() {
    return {
      imageList: [
        
      ],
      url2: "https://fuasdfss10.elemecdn.com/1/34/19aa98b1fcb2781c4fba33d850549jpeg.jpeg",
      url3: "https://fuss10.elemecdn.com/0/6f/e35ff375812e6b0020b6b4e8f9583jpeg.jpeg",
      url4: "https://fuss10.elemecdn.com/9/bb/e27858e973f5d7d3904835f46abbdjpeg.jpeg",
    };
  },
  computed: {
    scrollContainer() {
      return this.$refs.el_scrollbar?.wrap;
    },
  },
  methods: {
    downloadImage(url) {
      if (!url) return;

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
        // 添加跳转方法
    jumpToTutorial(section) {
      // 修改路径为根路径
      window.open(`/tutorial.html#${section}`, '_blank');
    },

    
    async checkUserInfo() {
      try {
        const response = await get('check_user_info')  // 使用 get 方法
        if (response.data.success) {  // 注意：response.data 而不是 response.json()
          ElMessage.success(`当前用户ID: ${response.data.user_id}`);
        } else {
          ElMessage.warning(response.data.message || '获取用户信息失败');
        }
      } catch (error) {
        ElMessage.error('请求失败：' + error.message);
      }
    }
  }
};
</script>

<style scoped>
.photo-gallery {
  padding: 20px;
  width: 100%;
  max-width: 1200px;
  margin: 0 auto;
}

.gallery-row {
  display: flex;
  justify-content: space-between;
  margin-bottom: 30px;
  gap: 20px;
}

.photo-item {
  flex: 1;
  text-align: center;
}

.photo-item h3 {
  margin-bottom: 10px;
  color: #333;
}

.el-image {
  width: 100%;
  height: 300px;
  border-radius: 8px;
  box-shadow: 0 2px 12px 0 rgba(0, 0, 0, 0.1);
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

.title-container {
  display: flex;
  align-items: center;
  justify-content: center;
  gap: 8px;
}

.title-container h3 {
  margin: 0;
}
</style>
