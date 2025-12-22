import { ElMessage } from 'element-plus';
import router from '../router';

export const autoLogoutManager = {
  timer: null,
  // timeoutDuration: 3 * 1000, // 3s
  timeoutDuration: 12 * 60 * 60 * 1000, // 12h
  boundResetTimer: null,
  isInitialized: false, // 添加标志位追踪初始化状态

  init() {
    if (!this.isInitialized) {
      this.boundResetTimer = () => this.resetTimer();
      
      // 添加事件监听
      document.addEventListener('mousemove', this.boundResetTimer);
      document.addEventListener('keypress', this.boundResetTimer);
      document.addEventListener('click', this.boundResetTimer);
      
      // 初始化计时器
      this.resetTimer();
      this.isInitialized = true;
    }
  },

  resetTimer() {
    // 只在已登录状态下重置计时器
    if (localStorage.getItem('user_id')) {
      if (this.timer) {
        clearTimeout(this.timer);
      }
      this.timer = setTimeout(() => {
        this.logout();
      }, this.timeoutDuration);
    } else {
      // 如果未登录，清理所有状态
      this.removeEventListeners();
    }
  },

  removeEventListeners() {
    if (this.boundResetTimer) {
      document.removeEventListener('mousemove', this.boundResetTimer);
      document.removeEventListener('keypress', this.boundResetTimer);
      document.removeEventListener('click', this.boundResetTimer);
    }
    if (this.timer) {
      clearTimeout(this.timer);
      this.timer = null;
    }
    this.isInitialized = false;
  },

  logout() {
    // 清除登录信息
    localStorage.removeItem('user_id');
    localStorage.removeItem('user_name');
    
    // 清理事件监听和计时器
    this.removeEventListeners();

    // 提示用户
    ElMessage({
      message: '由于长时间未操作，已自动登出',
      type: 'warning'
    });

    // 跳转到登录页
    if (router.currentRoute.value.path !== '/login') {
      router.push('/login');
    }
  }
};