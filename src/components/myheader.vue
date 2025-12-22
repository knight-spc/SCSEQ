<template>
  <div class="header-container">
    <el-page-header @back="goBack">
      <template #content>
        <span class="text-large font-500 mr-3"> {{ this.$route.meta.title }} </span>
      </template>
    </el-page-header>

    <div class="user-info">
      <el-dropdown @command="handleCommand">
        <span class="el-dropdown-link">
          <el-avatar :size="33" icon="el-icon-user-solid" :src="circleUrl" style="margin-right: 8px; vertical-align: middle;"></el-avatar>
          <span style="vertical-align: middle; font-size: 16px;">{{ username }}</span>
          <el-icon class="el-icon--right"><arrow-down /></el-icon>
        </span>
        <template #dropdown>
          <!-- 根据登录状态显示不同的下拉菜单 -->
          <el-dropdown-menu v-if="isLoggedIn">
            <el-dropdown-item command="profile" class="menu-item">
              <el-icon><User /></el-icon>个人中心
            </el-dropdown-item>
            <el-dropdown-item command="logout" class="menu-item">
              <el-icon><SwitchButton /></el-icon>退出登录
            </el-dropdown-item>
          </el-dropdown-menu>
          <el-dropdown-menu v-else>
            <div class="benefits-item">
              <div class="login-benefits">
                <p class="benefits-title">登录后可享受以下权益：</p>
                <ul class="benefits-list">
                  <li><el-icon><Files /></el-icon> 更全面的功能</li>
                  <li><el-icon><DataLine /></el-icon> 更丰富的可视化</li>
                  <li><el-icon><Connection /></el-icon> 更高效的创作环境</li>
                </ul>
              </div>
            </div>
            <el-divider class="benefits-divider" />
            <el-dropdown-item command="login" class="login-button">
              <el-icon><User /></el-icon>登录/注册
            </el-dropdown-item>
          </el-dropdown-menu>
        </template>
      </el-dropdown>
    </div>
  </div>
</template>

<script>
import { ArrowDown, User, SwitchButton, Document, ChatLineRound, Connection } from '@element-plus/icons-vue'
import defaultAvatar from '@/assets/images/user.png';

export default {
  components: {
    ArrowDown,
    User,
    SwitchButton,
    Document,
    ChatLineRound,
    Connection
  },
  data() {
    return {
      username: '未登录',
      circleUrl: defaultAvatar,
      isLoggedIn: false,
    };
  },
  mounted() {
    // 检查登录状态
    const userId = localStorage.getItem('user_id');
    this.isLoggedIn = !!userId;
    if (this.isLoggedIn) {
      this.username = localStorage.getItem('user_name');
    }
  },
  methods: {
    goBack() {
      this.$router.back();
      console.log("go back");
    },
    handleCommand(command) {
      if (command === 'logout') {
        this.logout();
      } else if (command === 'profile') {
        this.goToProfile();
      } else if (command === 'login') {
        this.$router.push('/login');
      }
    },
    logout() {
      // 清除本地存储的用户信息
      localStorage.removeItem('user_id');
      localStorage.removeItem('user_name');
      // 跳转到登录页
      this.$router.push('/login');
      console.log("User logged out");
    },
    goToProfile() {
      // 跳转到个人中心页面，假设路由为 /profile
      this.$router.push('/profile');
      console.log("Go to profile page");
    }
  },
};
</script>

<style scoped>
.header-container {
  display: flex;
  justify-content: space-between;
  align-items: center;
  padding: 10px 6% 10px 0px; /* 左边距20px，右边距5px */
  height: 60px;
  box-sizing: border-box;
  /* 添加负的上外边距，尝试将组件向上移动 */
  margin-top: -11px; /* 你可以调整这个负值的大小 */
}


.user-info {
  cursor: pointer;
}

.el-dropdown-link {
  display: flex;
  align-items: center;
  color: #409EFF; /* Element Plus 主题色 */
}

.el-dropdown-link:focus {
    outline: none; /* 移除焦点时的轮廓 */
}

/* 可以根据需要调整头像和用户名的样式 */
.el-avatar {
  background-color: #c0c4cc; /* 默认头像背景色 */
}

.benefits-item {
  cursor: default !important;
  padding: 12px 16px !important;
}

.benefits-item.el-dropdown-item.is-disabled {
  background-color: transparent;
  color: inherit;
}

.benefits-item:hover {
  background-color: transparent !important;
}

.login-benefits {
  min-width: 240px;
}

.benefits-title {
  margin: 0;
  padding: 5px 0 10px;
  color: #303133;
  font-weight: 600;
  font-size: 16px;
}

.benefits-list {
  list-style: none;
  padding: 0;
  margin: 0;

}

.benefits-list li {
  display: flex;
  align-items: center;
  padding: 8px 0;
  color: #606266;
  font-size: 14px;
  transition: all 0.3s;
}

.benefits-list li:hover {
  color: #409EFF;
}

.benefits-list .el-icon {
  margin-right: 10px;
  color: #409EFF;
  font-size: 18px;
}

.benefits-divider {
  margin: 4px 0 !important;
}

.login-button {
  text-align: center;
  color: #409EFF !important;
  font-weight: 500;
}

.login-button:hover {
  background-color: #ecf5ff !important;
}

.login-button .el-icon {
  margin-right: 6px;
  font-size: 15px;
}

.menu-item {
  display: flex;
  align-items: center;
  padding: 8px 16px;
  color: #606266;
  font-size: 14px;
  transition: all 0.3s;
  
}

.menu-item:hover {
  color: #409EFF !important;
  background-color: #ecf5ff !important;
}

.menu-item .el-icon {
  margin-right: 10px;
  color: #409EFF;
  font-size: 18px;
}
</style>
  