
<template>
  <div id="app">
    <router-view></router-view>
  </div>
</template>

<script>
import { autoLogoutManager } from '@/utils/autoLogout.js'

export default {
  name: 'App',
  mounted() {
    window.addEventListener('storage', this.handleStorageChange);
    
    if (localStorage.getItem('user_id')) {
      autoLogoutManager.init();
    }
  },
  methods: {
    handleStorageChange(event) {
      if (event.key === 'user_id' && event.newValue) {
        autoLogoutManager.init();
      }
      else if (event.key === 'user_id' && !event.newValue) {
        autoLogoutManager.removeEventListeners();
      }
    }
  },
  beforeUnmount() {
    window.removeEventListener('storage', this.handleStorageChange);
    autoLogoutManager.removeEventListeners();
  }
}
</script>

<style>
#app {
  width: 100%;
  height: 100%;
}
</style>
