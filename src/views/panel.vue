<!--仪表板-->

<template>
  <div id="app">
    <el-container>
      <el-aside width="200px">
        <sidebar></sidebar>
      </el-aside>
      <el-container >
        <el-header class="el-header">
          <myheader></myheader>
        </el-header>
        <el-divider style="margin: 0 0 1px" />
        <el-main style="padding-top: 5px;padding-bottom:0px;"><el-scrollbar ref="el_scrollbar" style=" height: 600px;overflow-y: auto; "><router-view></router-view>
          </el-scrollbar></el-main>
        <el-footer class="el-footer">
          <myfooter></myfooter>
        </el-footer>
      </el-container>
    </el-container>
  </div>
</template>
  
<script>
import sidebar from "../components/sidebar.vue";
import myheader from "../components/myheader.vue";
import myfooter from "../components/myfooter.vue";
import SanJiKaiShu from '../assets/font/font.css';

import { useRoute } from "vue-router";

export default {
  components: {
    sidebar,
    myheader,
    myfooter,
  },
  data() {
    return {
      show_sidebar: false,
    };
  },
  methods: {},
  mounted() { },
};
</script>
  
<style>
#app {
  font-family: Avenir, Helvetica, Arial, sans-serif;
  /* font-family: SanJiKaiShu; */
  font-weight: bold;
}

.el-header {
  --el-header-padding: 20px 20px 0 !important;
  /* background-color:; */
}
</style>
  
  
  
  