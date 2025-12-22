<template>
      <el-container>
        <el-header>
          <h2>注册页面</h2>
        </el-header>
        <el-main>
          <el-form :model="registrationForm" :rules="rules" ref="registrationForm" label-width="120px">
            <el-form-item label="用户名" prop="username">
              <el-input v-model="registrationForm.username"></el-input>
            </el-form-item>
            <el-form-item label="邮箱" prop="email">
              <el-input v-model="registrationForm.email"></el-input>
            </el-form-item>
            <el-form-item label="密码" prop="password">
              <el-input type="password" v-model="registrationForm.password"></el-input>
            </el-form-item>
            <el-form-item label="确认密码" prop="confirmPassword">
              <el-input type="password" v-model="registrationForm.confirmPassword"></el-input>
            </el-form-item>
            <el-form-item>
              <el-button type="primary" @click="submitForm('registrationForm')">注册</el-button>
              <el-button @click="resetForm('registrationForm')">重置</el-button>
            </el-form-item>
          </el-form>
        </el-main>
      </el-container>
  </template>
  
  <script>
  export default {
    data() {
      return {
        registrationForm: {
          username: '',
          email: '',
          password: '',
          confirmPassword: ''
        },
        rules: {
          username: [
            { required: true, message: '请输入用户名', trigger: 'blur' },
            { min: 3, max: 20, message: '用户名长度在 3 到 20 个字符', trigger: 'blur' }
          ],
          email: [
            { required: true, message: '请输入邮箱地址', trigger: 'blur' },
            { type: 'email', message: '请输入正确的邮箱地址', trigger: ['blur', 'change'] }
          ],
          password: [
            { required: true, message: '请输入密码', trigger: 'blur' },
            { min: 6, message: '密码长度不能少于 6 个字符', trigger: 'blur' }
          ],
          confirmPassword: [
            { required: true, message: '请确认密码', trigger: 'blur' },
            { validator: this.validateConfirmPassword, trigger: 'blur' }
          ]
        }
      };
    },
    methods: {
      validateConfirmPassword(rule, value, callback) {
        if (value === '') {
          callback(new Error('请确认密码'));
        } else if (value !== this.registrationForm.password) {
          callback(new Error('两次输入密码不一致!'));
        } else {
          callback();
        }
      },
      submitForm(formName) {
        this.$refs[formName].validate((valid) => {
          if (valid) {
            alert('注册成功!');
          } else {
            console.log('error submit!!');
            return false;
          }
        });
      },
      resetForm(formName) {
        this.$refs[formName].resetFields();
      }
    }
  };
  </script>
  
  <style>
  #app {
    font-family: Avenir, Helvetica, Arial, sans-serif;
    -webkit-font-smoothing: antialiased;
    -moz-osx-font-smoothing: grayscale;
    text-align: center;
    color: #2c3e50;
    margin-top: 60px;
  }
  
  el-container {
    width: 50%;
    margin: 0 auto;
  }
  </style>