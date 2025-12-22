import { createApp } from 'vue'

import content from '../content.js';
import App from './App.vue'
import router from './router'
import store from './store'
import ECharts from 'vue-echarts'

import ElementPlus from 'element-plus'
import 'element-plus/dist/index.css'
import axios from 'axios'
import VueAxios from "vue-axios";
import naive from "naive-ui";

const app = createApp(App)
 
app.use(store).use(router).use(ElementPlus)
app.use(VueAxios, axios).use(router) 
app.use(naive);
app.provide('axios', app.config.globalProperties.axios) 
app.mount('#app')

app.config.globalProperties.$myContent = content ;

import * as ElementPlusIconsVue from '@element-plus/icons-vue'

for (const [key, component] of Object.entries(ElementPlusIconsVue)) {
  app.component(key, component)
}

const debounce = (fn, delay) => {
  let timer = null;
  return function () {
    let context = this;
    let args = arguments;
    clearTimeout(timer);
    timer = setTimeout(function () {
      fn.apply(context, args);
    }, delay);
  };
};

const _ResizeObserver = window.ResizeObserver;
window.ResizeObserver = class ResizeObserver extends _ResizeObserver {
  constructor(callback) {
    callback = debounce(callback, 16);
    super(callback);
  }
};

