<template>
  <el-card>
  <div class="form-container">
    <el-form
      label-position="left"
      label-width="100px"
      :model="form"
      style="max-width: 600px"
    >
      <el-form-item label="选择导出变量">
        <el-radio-group v-model="data_name">
          <el-radio label="temp">温度</el-radio>
          <el-radio label="light">光照</el-radio>
          <el-radio label="humidity">湿度</el-radio>
          <el-radio label="CO2">CO2</el-radio>
          <el-radio label="weight">湿度</el-radio>
        </el-radio-group>
      </el-form-item>
      <el-form-item label="起始时间">
        <el-date-picker
          v-model="start_time"
          type="datetime"
          placeholder="选择一个时间"
          format="YYYY/MM/DDTHH:mm:ss"
          value-format="YYYY-MM-DDTHH:mm:ss"
        />
      </el-form-item>
      <el-form-item label="结束时间">
        <el-date-picker
          v-model="end_time"
          type="datetime"
          placeholder="选择一个时间"
          format="YYYY/MM/DDTHH:mm:ss"
          value-format="YYYY-MM-DDTHH:mm:ss"
        />
      </el-form-item>
      <el-form-item label="数据条数">
        <el-input-number
          v-model="limit"
          :min="1"
          :max="1000"
          id="limit"
        ></el-input-number>
      </el-form-item>
      <el-form-item>
        <el-button type="primary" @click="onSubmit" :loading="isLoading">导出数据</el-button>
      </el-form-item>
    </el-form>
  </div>
</el-card>
</template>
  
  <script>
import { reactive } from "vue";
// import { generate_excel } from "../request/api.js";
import { ElMessage } from "element-plus";

export default {
  data() {
    return {
      data_name: "temp",
      start_time: "",
      end_time: "",
      limit: 100,
      isLoading:false
    };
  },
  methods: {
    onSubmit() {
      this.loadingButton();
      generate_excel(this.data_name, this.start_time, this.end_time, this.limit)
        .then((res) => {
          const blob = new Blob([res.data], {
            type: "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
          });
          const url = URL.createObjectURL(blob);
          ElMessage({
            message: "导出成功",
            type: "success",
          });
          window.open(url, "_blank");
        })
        .catch((err) => {
          console.log(err);
          ElMessage.error(JSON.stringify(err));
        });
    },
    loadingButton(){
      this.isLoading = true;
      const timer = setTimeout(() => {
          clearTimeout(timer);
          this.isLoading = false;
        }, 500);
    }
  },
};
</script>
  
  <style scoped>
.form-container {
  display: flex;
  justify-content: center;
  align-items: center;
}
</style>
  