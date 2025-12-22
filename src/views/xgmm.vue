<template>
  <el-card>
    <div class="password-container">
      <el-form :model="form" style="max-width: 600px;">
        <el-form-item label="当前密码" prop="currentPassword">
          <el-input v-model="form.currentPassword" show-password></el-input>
        </el-form-item>
        <el-form-item label="新密码" prop="newPassword">
          <el-input v-model="form.newPassword" show-password></el-input>
        </el-form-item>
        <el-form-item label="确认新密码" prop="confirmNewPassword">
          <el-input v-model="form.confirmNewPassword" show-password></el-input>
        </el-form-item>
        <el-form-item>
          <el-button type="primary" @click="onSubmit">提交</el-button>
          <el-button @click="onCancel">取消</el-button>
        </el-form-item>
      </el-form>
    </div>
  </el-card>
  </template>
  
  <script>
  import { reactive } from 'vue';
  import { ElForm, ElFormItem, ElInput, ElButton, ElMessage } from 'element-plus';
  
  export default {
    components: {
      ElForm,
      ElFormItem,
      ElInput,
      ElButton,
    },
    data() {
      return {
        form: reactive({
          currentPassword: '',
          newPassword: '',
          confirmNewPassword: '',
        }),
      };
    },
    methods: {
      onSubmit() {
        if (this.form.newPassword !== this.form.confirmNewPassword) {
          ElMessage.error('两次输入不匹配！');
          return;
        }
        
        // Call API to update password
        
        ElMessage.success('修改成功');
      },
      onCancel() {
        // Navigate back to previous page
      },
    },
  };
  </script>
  
  <style scoped>
  .password-container {
    display: flex;
    justify-content: center;
    align-items: center;
  }
  </style>
  