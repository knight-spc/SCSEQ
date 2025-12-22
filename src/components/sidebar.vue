<!--该文件表示右侧菜单栏-->

<template>
  <el-col>
   <div>
       <a href="/"><img style="width:90%" src="../assets/images/scseq_logo.png"></a>
   </div>
 <el-menu
   class="el-menu-vertical-demo"
   :default-active="this.$route.path"
   :default-openeds="['1', '2']"
   @open="handleOpen"
   @close="handleClose"
   @select="handleMenuSelect"
 >
   <el-sub-menu index="1">
     <template #title>
       <el-icon><User /></el-icon>
       <span>系统用户</span>
     </template>
       <el-menu-item index="yhgl">用户管理</el-menu-item>
       <!-- <el-menu-item index="yhxx">用户信息</el-menu-item> -->
   </el-sub-menu>
   <el-sub-menu index="2">
       <template #title>
         <el-icon><Stopwatch /></el-icon>
         <span>测序分析</span>
       </template>
     <el-menu-item index="upload">项目管理</el-menu-item>
     <el-menu-item index="sdcs">设定参数</el-menu-item>
     <el-menu-item index="jsjg">查看可视化结果</el-menu-item>
     <el-menu-item index="znjy">更多信息</el-menu-item>
     
   </el-sub-menu>
   <!-- <el-menu-item index="sjwh">
     <el-icon><DataAnalysis /></el-icon>
     <span>数据维护</span>
   </el-menu-item> -->
   <!-- <el-menu-item index="xgmm">
     <el-icon><VideoCamera /></el-icon>
     <span>修改密码</span>
   </el-menu-item>
   <el-menu-item index="znfx">
     <el-icon><Odometer /></el-icon>
     <span>智能分析</span>
   </el-menu-item> -->
   <!-- <el-menu-item index="tc">
     <el-icon><Back /></el-icon>
     <span>退出</span>
   </el-menu-item> -->
   
 </el-menu>
</el-col>
</template>

<script>
import {
User,
Stopwatch,
DataAnalysis,
VideoCamera,
Back,
Odometer
} from "@element-plus/icons-vue";

export default {
components:{
   User,
   Stopwatch,
   DataAnalysis,
   VideoCamera,
   Back,
   Odometer
},
name:'sidebar',
methods: {
handleMenuSelect(index) {
 this.$router.push(index)
},
handleOpen(key, keyPath) {
 console.log(key, keyPath);
},
handleClose(key, keyPath) {
 console.log(key, keyPath);
},
},
};
</script>

<style scoped>

</style>