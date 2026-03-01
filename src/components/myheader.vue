<template>
  <div class="header-container">
    <el-page-header @back="goBack">
      <template #content>
        <span class="text-large font-500 mr-3">
          {{ $t($route.meta.title) }}
        </span>
      </template>
    </el-page-header>

    <!-- 右侧区域 -->
    <div class="right-actions">
      <!-- 语言切换 -->
      <el-dropdown @command="changeLanguage">
        <span class="lang-switch">
          <el-icon><Connection /></el-icon>
          <span class="lang-text">
            {{ currentLangLabel }}
          </span>
        </span>

        <template #dropdown>
          <el-dropdown-menu>
            <el-dropdown-item command="zh">中文</el-dropdown-item>
            <el-dropdown-item command="en">English</el-dropdown-item>
          </el-dropdown-menu>
        </template>
      </el-dropdown>

      <!-- 用户信息 -->
      <div class="user-info">
        <el-dropdown @command="handleCommand">
          <span class="el-dropdown-link">
            <el-avatar :size="33" :src="circleUrl" style="margin-right: 8px" />
            <span style="font-size: 16px">{{ username }}</span>
            <el-icon class="el-icon--right"><ArrowDown /></el-icon>
          </span>
          <template #dropdown>
            <!-- 根据登录状态显示不同的下拉菜单 -->
            <el-dropdown-menu v-if="isLoggedIn">
              <el-dropdown-item command="profile" class="menu-item">
                <el-icon><User /></el-icon>
                {{ $t("user.profile") }}
              </el-dropdown-item>

              <el-dropdown-item command="logout" class="menu-item">
                <el-icon><SwitchButton /></el-icon>
                {{ $t("user.logout") }}
              </el-dropdown-item>
            </el-dropdown-menu>

            <el-dropdown-menu v-else>
              <div class="benefits-item">
                <div class="login-benefits">
                  <p class="benefits-title">
                    {{ $t("loginBenefits.title") }}
                  </p>

                  <ul class="benefits-list">
                    <li>
                      <el-icon><Files /></el-icon>
                      {{ $t("loginBenefits.fullFeatures") }}
                    </li>
                    <li>
                      <el-icon><DataLine /></el-icon>
                      {{ $t("loginBenefits.richVisualization") }}
                    </li>
                    <li>
                      <el-icon><Connection /></el-icon>
                      {{ $t("loginBenefits.efficientWorkflow") }}
                    </li>
                  </ul>
                </div>
              </div>

              <el-divider class="benefits-divider" />

              <el-dropdown-item command="login" class="login-button">
                <el-icon><User /></el-icon>
                {{ $t("login.loginOrRegister") }}
              </el-dropdown-item>
            </el-dropdown-menu>
          </template>
        </el-dropdown>
      </div>
    </div>
  </div>
</template>

<script>
import {
  ArrowDown,
  User,
  SwitchButton,
  Document,
  ChatLineRound,
  Connection,
} from "@element-plus/icons-vue";
import defaultAvatar from "@/assets/images/user.png";

export default {
  components: {
    ArrowDown,
    User,
    SwitchButton,
    Document,
    ChatLineRound,
    Connection,
  },
  data() {
    return {
      username: "未登录",
      circleUrl: defaultAvatar,
      isLoggedIn: false,
      currentLang: this.$i18n.locale || "zh",
    };
  },
  computed: {
    currentLangLabel() {
      return this.currentLang === "zh" ? "中文" : "EN";
    },
  },

  mounted() {
    const userId = localStorage.getItem("user_id");
    this.isLoggedIn = !!userId;
    if (this.isLoggedIn) {
      this.username = localStorage.getItem("user_name");
    }

    // 初始化语言
    const savedLang = localStorage.getItem("lang");
    if (savedLang) {
      this.currentLang = savedLang;
      this.$i18n.locale = savedLang;
    }
  },

  methods: {
    goBack() {
      this.$router.back();
      console.log("go back");
    },
    changeLanguage(lang) {
      this.currentLang = lang;
      this.$i18n.locale = lang;
      localStorage.setItem("lang", lang);
    },
    handleCommand(command) {
      if (command === "logout") {
        this.logout();
      } else if (command === "profile") {
        this.goToProfile();
      } else if (command === "login") {
        this.$router.push("/login");
      }
    },
    logout() {
      // 清除本地存储的用户信息
      localStorage.removeItem("user_id");
      localStorage.removeItem("user_name");
      // 跳转到登录页
      this.$router.push("/login");
      console.log("User logged out");
    },
    goToProfile() {
      // 跳转到个人中心页面，假设路由为 /profile
      this.$router.push("/profile");
      console.log("Go to profile page");
    },
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
  color: #409eff; /* Element Plus 主题色 */
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
  color: #409eff;
}

.benefits-list .el-icon {
  margin-right: 10px;
  color: #409eff;
  font-size: 18px;
}

.benefits-divider {
  margin: 4px 0 !important;
}

.login-button {
  text-align: center;
  color: #409eff !important;
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
  color: #409eff !important;
  background-color: #ecf5ff !important;
}

.menu-item .el-icon {
  margin-right: 10px;
  color: #409eff;
  font-size: 18px;
}
.right-actions {
  display: flex;
  align-items: center;
  gap: 20px;
}

.lang-switch {
  display: flex;
  align-items: center;
  cursor: pointer;
  color: #606266;
  font-size: 14px;
}

.lang-switch:hover {
  color: #409eff;
}

.lang-text {
  margin-left: 6px;
}
</style>
