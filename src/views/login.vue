<template>
  <div class="login-page">
    <img src="../assets/images/scseq2.png" class="background-image" />
    <div class="login-box">
      <h1>单细胞转录组测序分析云平台</h1>
      <el-form v-if="!isRegister" :model="loginForm" :rules="rules" ref="loginForm">
        <el-form-item label="用户名" prop="username">
          <el-input v-model="loginForm.username" placeholder="请输入用户名"></el-input>
        </el-form-item>
        <el-form-item label="密码" prop="password">
          <el-input type="password" v-model="loginForm.password" placeholder="请输入密码"></el-input>
        </el-form-item>
        <el-form-item>
          <el-button type="primary" @click="submitForm" :loading="loading" style="width: 100%">登录</el-button>
        </el-form-item>
        <div class="form-footer">
          <span>首次登陆？</span>
          <el-button type="text" @click="toggleForm">去注册</el-button>
        </div>
      </el-form>

      <el-form v-else :model="registerForm" :rules="registerRules" ref="registerForm">
        <el-form-item label="用户名" prop="username">
          <el-input v-model="registerForm.username" placeholder="请输入用户名"></el-input>
        </el-form-item>
        <el-form-item label="密码" prop="password">
          <el-input type="password" v-model="registerForm.password" placeholder="请输入密码"></el-input>
        </el-form-item>
        <el-form-item label="手机号" prop="phone">
          <el-input v-model="registerForm.phone" placeholder="请输入手机号"></el-input>
        </el-form-item>
        <el-form-item>
          <el-button type="primary" @click="submitRegister" :loading="loading" style="width: 100%">注册</el-button>
        </el-form-item>
        <div class="form-footer">
          <span>已有账号？</span>
          <el-button type="text" @click="toggleForm">去登录</el-button>
        </div>
      </el-form>
    </div>
  </div>
</template>

<script>
import { ElMessage } from "element-plus";
import axios from 'axios';

export default {
  data() {
    return {
      isRegister: false,
      loading: false,
      loginForm: {
        username: "",
        password: "",
      },
      registerForm: {
        username: "",
        password: "",
        phone: ""
      },
      rules: {
        username: [
          { required: true, message: '请输入用户名', trigger: 'blur' },
          { pattern: /^[a-zA-Z0-9_]+$/, message: '用户名只能包含字母、数字和下划线' }
        ],
        password: [
          { required: true, message: '请输入密码', trigger: 'blur' },
          { pattern: /^[a-zA-Z0-9!@#$%^&*()_+]+$/, message: '密码只能包含字母、数字和特殊字符' }
        ]
      },
      registerRules: {
        username: [
          { required: true, message: '请输入用户名', trigger: 'blur' },
          { pattern: /^[a-zA-Z0-9_]+$/, message: '用户名只能包含字母、数字和下划线' }
        ],
        password: [
          { required: true, message: '请输入密码', trigger: 'blur' },
          { pattern: /^[a-zA-Z0-9!@#$%^&*()_+]+$/, message: '密码只能包含字母、数字和特殊字符' }
        ],
        phone: [
          { required: true, message: '请输入手机号', trigger: 'blur' },
          { pattern: /^1[3-9]\d{9}$/, message: '请输入正确的手机号' }
        ]
      },
      // 添加自动登出时间设置（30分钟）
      autoLogoutTime: 3 * 1000,
    };
  },
  methods: {
    toggleForm() {
      this.isRegister = !this.isRegister;
    },
    submitForm() {
      this.$refs.loginForm.validate(async (valid) => {
        if (valid) {
          this.loading = true;
          try {
            const formData = new FormData();
            formData.append('user_name', this.loginForm.username);
            formData.append('password', this.loginForm.password);
            
            const response = await axios.post('/api/login', formData, {
              headers: {
                'Content-Type': 'application/x-www-form-urlencoded'
              }
            });
            
            if (response.data.status === 'success') {
              localStorage.setItem('user_id', response.data.user_id);
              localStorage.setItem('user_name', response.data.user_name);
              localStorage.setItem('login_time', Date.now());
              
              ElMessage({
                message: '登录成功',
                type: 'success'
              });
              
              this.$router.push('/');
            }
          } catch (error) {
            let errorMessage = '登录失败';
            if (error.response) {
              console.error('错误响应:', error.response);
              errorMessage = error.response.data?.message || '服务器错误，请稍后重试';
            } else if (error.request) {
              console.error('请求错误:', error.request);
              errorMessage = '网络连接失败，请检查网络设置';
            } else {
              console.error('其他错误:', error.message);
              errorMessage = '发生未知错误，请重试';
            }
            ElMessage.error(errorMessage);
          } finally {
            this.loading = false;
          }
        }
      });
    },
    
    async submitRegister() {
      this.$refs.registerForm.validate(async (valid) => {
        if (valid) {
          this.loading = true;
          const formData = new FormData();
          formData.append('user_name', this.registerForm.username);
          formData.append('password', this.registerForm.password);
          formData.append('phone', this.registerForm.phone);

          try {
            const response = await axios.post('/api/register', formData);
            
            if (response.data.status === 'success') {
              ElMessage({
                message: '注册成功',
                type: 'success'
              });
              this.isRegister = false; // 切换到登录页
            }
          } catch (error) {
            let errorMessage = '注册失败';
            if (error.response && error.response.data) {
              errorMessage = error.response.data.message;
            }
            ElMessage.error(errorMessage);
          } finally {
            this.loading = false;
          }
        }
      });
    },
    

  },
  
};
</script>

<style>
.login-page {
  height: 100vh;
  width: 100vw;
  overflow: hidden;
  position: relative;
  display: flex;
  justify-content: center;
  align-items: center;
  margin: 0;
  padding: 0;
}

.background-image {
  position: absolute;
  top: 0;
  left: 0;
  width: 100%;
  height: 100%;
  object-fit: cover;
  z-index: -1;
}

.login-box {
  width: 400px;
  padding: 40px;
  background: rgba(255, 255, 255, 0.95);
  border-radius: 8px;
  box-shadow: 0 8px 16px rgba(0, 0, 0, 0.1);
  backdrop-filter: blur(5px);
  position: relative;
  z-index: 1;
}

.login-box h1 {
  color: #333;
  font-size: 24px;
  text-align: center;
  margin-bottom: 30px;
  font-weight: 500;
}

.form-footer {
  text-align: center;
  margin-top: 20px;
  color: #606266;
}

.el-form-item {
  margin-bottom: 25px;
}

.el-button--text {
  color: #409EFF;
  margin-left: 8px;
  font-weight: 500;
}

.el-input {
  height: 40px;
}

/* 添加一些动画效果 */
.el-button {
  transition: all 0.3s ease;
}

.el-button:hover {
  transform: translateY(-2px);
  box-shadow: 0 2px 8px rgba(0, 0, 0, 0.1);
}

/* 确保页面占满整个视窗 */
body {
  margin: 0;
  padding: 0;
  overflow: hidden;
}
</style>
  