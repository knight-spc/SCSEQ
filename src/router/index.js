import { createRouter, createWebHistory } from 'vue-router'
import panel from "@/views/panel.vue";
import Welcome from '../views/welcome.vue'; // 确保路径正确

// ()=>import('../views/panel.vue')
const routes = [
  {
    path: '/',
    name: 'home',
    component: panel,
    children:[
      {
        path: '',
        name: 'welcome',
        component: () => import(/* webpackChunkName: "about" */ '../views/welcome.vue'),
        meta:{ title:"欢迎页" }
      },
      {
        path: '/yhgl',
        name: 'yhgl',
        component: () => import(/* webpackChunkName: "about" */ '../views/yhgl.vue'),
        meta:{ title:"用户管理" }
      },
      {
        path: '/xgmm',
        name: 'xgmm',
        component: () => import(/* webpackChunkName: "about" */ '../views/xgmm.vue'),
        meta:{ title:"修改密码" }
      },
      
      {
        path: '/upload',
        name: 'upload',
        component: () => import(/* webpackChunkName: "about" */ '../views/upload.vue'),
        meta:{ title:"项目管理" }
      },
      
     
      {
        path: '/sjwh',
        name: 'sjwh',
        component: () => import(/* webpackChunkName: "about" */ '../views/sjwh.vue'),
        meta:{ title:"数据维护" }
      },
      
      {
        path: '/sdcs',
        name: 'sdcs',
        component: () => import(/* webpackChunkName: "about" */ '../views/sdcs.vue'),
        meta:{ title:"设定参数" }
      },
      {
        path: '/sysxc',
        name: 'sysxc',
        component: () => import(/* webpackChunkName: "about" */ '../views/sysxc.vue'),
        meta:{ title:"实验室现场" }
      },
      {
        path: '/znjy',
        name: 'znjy',
        component: () => import(/* webpackChunkName: "about" */ '../views/znjy.vue'),
        meta:{ title:"智能建议" }
      },
      {
        path: '/jsjg',
        name: 'jsjg',
        component: () => import(/* webpackChunkName: "about" */ '../views/jsjg.vue'),
        meta:{ title:"计算结果" }
      },
      {
        path: '/znfx',
        name: 'znfx',
        component: () => import(/* webpackChunkName: "about" */ '../views/znfx.vue'),
        meta:{ title:"智能分析" }
      },
    ]
  },
  {
    path: '/login',
    name: 'login',
    component: () => import(/* webpackChunkName: "about" */ '../views/login.vue'),
    // meta: { hideSidebar: true }
  },
  {
    path: '/register',
    name: 'register',
    component: () => import(/* webpackChunkName: "about" */ '../views/register.vue'),
    meta: { hideSidebar: true }
  },
  
]
import { autoLogoutManager } from '@/utils/autoLogout.js'
const router = createRouter({
  history: createWebHistory(process.env.BASE_URL),
  routes
})
// import { createRouter, createWebHashHistory } from 'vue-router';
// const router = createRouter({
//   history: createWebHashHistory(),
//   routes
// });
// 添加全局前置守卫
router.beforeEach((to, from, next) => {
  const isLoggedIn = !!localStorage.getItem('user_id')
  
  if (to.path !== '/login' && !isLoggedIn) {
    next({
      path: '/login',
      query: { redirect: to.fullPath }
    })
  } else {
    // 如果用户已登录，确保初始化自动登出功能
    if (isLoggedIn) {
      autoLogoutManager.init()
    }
    next()
  }
})

export default router

