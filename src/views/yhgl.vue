<template>
  <el-table :data="tableData" stripe style="width: 100%">
    <el-table-column type="index" width="50" label="No." />
    <el-table-column prop="name" label="账户名称" />
    <el-table-column prop="state" label="备注" />
    <el-table-column fixed="right" label="操作">
      <template #default="{ row }">
        <el-button type="danger" size="small" @click="handleDelete(row)"
          >删除</el-button
        >
        <el-button type="primary" size="small" @click="handleEdit(row)"
          >编辑</el-button
        >
      </template>
    </el-table-column>
  </el-table>

  <div class="add">
    <el-button type="primary" @click="handleAdd">增加</el-button>
  </div>

  <el-dialog
    v-model="dialogVisible"
    title="用户详情"
    :visible="dialogVisible"
    @close="resetForm"
  >
    <el-form :model="formData" ref="form" label-width="80px">
      <el-form-item label="账户名称">
        <el-input v-model="formData.name"></el-input>
      </el-form-item>
      <el-form-item label="备注">
        <el-input v-model="formData.state"></el-input>
      </el-form-item>
    </el-form>
    <div>
      <el-button @click="dialogVisible = false">取 消</el-button>
      <el-button type="primary" @click="save">确 定</el-button>
    </div>
  </el-dialog>
</template>

<script>
import { ElMessage } from "element-plus";

export default {
  data() {
    return {
      tableData: [
        {
          id: 1,
          name: "老贾",
          state: "哈肥",
          city: "Los Angeles",
          address: "No. 189, Grove St, Los Angeles",
          zip: "CA 90036",
          tag: "Home",
        },
        {
          id: 2,
          name: "小邵",
          state: "上海",
          city: "Los Angeles",
          address: "No. 189, Grove St, Los Angeles",
          zip: "CA 90036",
          tag: "Office",
        },
        {
          id: 3,
          name: "小孙",
          state: "安徽",
          city: "Los Angeles",
          address: "No. 189, Grove St, Los Angeles",
          zip: "CA 90036",
          tag: "Home",
        },
        {
          id: 4,
          name: "小王",
          state: "山东",
          city: "Los Angeles",
          address: "No. 189, Grove St, Los Angeles",
          zip: "CA 90036",
          tag: "Office",
        },
      ],
      dialogVisible: false,
      formData: {},
    };
  },
  methods: {
    handleDelete(row) {
      const index = this.tableData.indexOf(row);
      this.tableData.splice(index, 1);
      ElMessage({
        message: "删除成功",
        type: "success",
      });
    },
    handleEdit(row) {
      this.dialogVisible = true;
      this.formData = { ...row };
    },
    save() {
      if (this.formData.id) {
        const index = this.tableData.findIndex(
          (item) => item.id === this.formData.id
        );
        this.tableData.splice(index, 1, this.formData);
      } else {
        const newId = this.tableData.length + 1;
        this.tableData.push({ ...this.formData, id: newId });
      }
      this.resetForm();
      this.dialogVisible = false;
      ElMessage({
        message: "成功",
        type: "success",
      });
    },
    resetForm() {
      this.formData = {};
      this.$refs.form.resetFields();
    },
    handleAdd() {
      this.dialogVisible = true;
    },
  },
};
</script>

<style>
.add {
  display: flex;
  justify-content: center;
  margin-top: 20px;
}
</style>
