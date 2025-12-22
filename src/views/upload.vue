<template>
  <el-card class="rounded-card">
    <div class="projects-header">
      <h2>我的项目</h2>
      <div class="button-group">
        <el-button plain @click="dialogFormVisible = true">新建项目</el-button>
        <el-tooltip
          class="box-item"
          effect="dark"
          content="介绍手册"
          placement="top"
        >
          <el-button type="" text size="small" class="help-button" @click="jumpToTutorial('data-upload')">
            <el-icon><QuestionFilled /></el-icon>
          </el-button>
        </el-tooltip>
      </div>
    </div>


    <!-- 项目列表 -->
    <div class="projects-grid">
      <el-card 
        v-for="item in items" 
        :key="item.id" 
        class="project-card" 
        shadow="hover"
        @click="handleCardClick(item)"
      >
        <template #header>
          <div class="project-header">
            <div class="project-info">
              <span class="project-name">{{ item.name }}</span>
              <span class="project-time">{{ item.created_time }}</span>
            </div>
            <el-tag size="small" :type="getSpeciesTagType(item.species)">
              {{ item.species }}
            </el-tag>
          </div>
        </template>
        <div class="project-content">
          <p class="project-desc">{{ item.notes || "暂无描述" }}</p>
          <div class="project-actions">
            <el-tooltip
              class="box-item"
              effect="dark"
              content="查看"
              placement="top"
            >
              <el-button type="primary" text class="action-button" @click.stop="viewProject(item)">
                <el-icon><View /></el-icon>
              </el-button>
            </el-tooltip>
            <el-tooltip
              class="box-item"
              effect="dark"
              content="删除"
              placement="top"
            >
              <el-button type="danger" text class="action-button" @click.stop="deleteProject(item)">
                <el-icon><Delete /></el-icon>
              </el-button>
            </el-tooltip>
          </div>
        </div>
      </el-card>
    </div>

    <!-- 原有的dialog保持不变 -->
    <el-dialog v-model="dialogFormVisible" title="上传数据" width="80%" height="100%">
      <el-row :gutter="30">
        <el-col :span="12">
          <el-form ref="form" :model="form" label-width="100px">
            <el-form-item label="项目名称：">
              <el-input v-model="form.name"></el-input>
            </el-form-item>

            <el-form-item label="物种信息：">
              <el-select v-model="form.species" placeholder="请选择物种信息">
                <el-option label="human" value="human"></el-option>
                <el-option label="mouse" value="mouse"></el-option>
                <el-option label="其他" value="other_species"></el-option>
              </el-select>
            </el-form-item>

            <el-form-item label="项目备注：">
              <el-input type="textarea" v-model="form.desc"></el-input>
            </el-form-item>
          </el-form>

          <div class="block">
            <div class="insert"></div>
            <div class="container">
              <div class="panel">
                <div class="table_head">请在这里上传测序数据</div>
                <el-upload
                  class="upload-demo"
                  drag
                  action="/api/upload"
                  name="img"
                  id="123"
                  multiple
                  accept=".zip"
                  :on-success="handleSuccess"
                  :on-error="handleError"
                  style="flex-direction: column; display: flex; justify-content: center"
                >
                  <el-icon class="el-icon--upload"><upload-filled /></el-icon>
                  <div class="el-upload__text">
                    Drop file here or <em>click to upload</em>
                  </div>
                  <template #tip>
                    <div class="el-upload__tip" style="text-align: center">
                      注意：打包前请保证文件格式正确
                    </div>
                  </template>
                </el-upload>
              </div>
            </div>
            <div class="insert"></div>
          </div>
        </el-col>
        <el-col :span="12">
          <el-row> 上传文件的格式要求： </el-row>
          <el-row>
            <div class="grid-content">
              <el-image :src="upload_example"> </el-image>
            </div>
          </el-row>
        </el-col>
      </el-row>

      <div slot="footer" class="dialog-footer">
        <el-button @click="dialogFormVisible = false">取 消</el-button>
        <el-button type="primary" @click="submitForm('form')">确 定</el-button>
      </div>
    </el-dialog>
    <!-- 任务信息弹窗 -->
    <el-dialog v-model="taskDialogVisible" title="任务信息" width="80%" :modal="true" :close-on-click-modal="false" style="max-width:1200px;">
      <!-- 树状图 -->
      <div ref="treeChart" v-loading="taskLoading" style="width: 100%; height: 600px; background-color: #ffffff; border-radius: 4px;"></div>
      
      <!-- 任务管理表格 -->
      <el-table :data="taskList" v-loading="taskLoading" style="width: 100%; margin-top: 20px;">
        <el-table-column prop="type" label="任务类型" width="130"/>
        <el-table-column prop="parameter" label="参数" />
        <el-table-column prop="status" label="运行状态" width="120">
          <template #default="scope">
            <el-tag
              v-if="scope.row.status === '运行中'"
              type="warning"
              effect="plain"
            >
              进行中
              <i class="el-icon-loading" style="margin-left:4px;"></i>
            </el-tag>
            <el-tag
              v-else
              type="success"
              effect="plain"
            >
              已完成
            </el-tag>
          </template>
        </el-table-column>
        <el-table-column prop="time" label="提交时间" width="180"/>
        <el-table-column label="操作" width="140">
          <template #default="scope">
            <el-button
              type="primary"
              size="small"
              :disabled="scope.row.status === '运行中'"
              style="margin-right:1px"
              @click="handleTaskView(scope.row)"
            >查看</el-button>
            <el-button
              type="danger"
              size="small"
              :disabled="scope.row.status === '运行中'"
              @click="handleTaskDelete(scope.row)"
            >删除</el-button>
          </template>
        </el-table-column>
      </el-table>
    </el-dialog>

  </el-card>
</template>
<script>
import * as echarts from 'echarts';
import { get, post } from "@/request/api";
import { UploadFilled, View, Delete } from '@element-plus/icons-vue';
export default {
  components: {
    UploadFilled,
    View,
    Delete,
  },
  data() {
    return {
      dialogFormVisible: false,
      upload_example: require("../assets/images/upload_example.png"),
      form: {
        name: "",
        species: "",
        desc: "",
        save_path: "",
        figure_path: "",
      },
      items: [], // 添加项目列表数据
      taskDialogVisible: false,
      taskList: [],
      taskLoading: false,
      treeData: null, // 新增
    };
  },

  created() {
    this.showItems();
  },

  methods: {
    // 上传数据
    handleSuccess(response, file) {
      this.num = response.task_id;
      this.form.save_path = response.save_path;
      this.form.figure_path = response.figure_path;

      this.$message({
        showClose: true,
        message: "数据上传成功",
        type: "success",
      });
    },
    handleError(response, file) {
      this.$message({
        showClose: true,
        message: "数据格式错误，请上传标准格式",
        type: "error",
      });
    },


    submitForm(formName) {
      const userId = localStorage.getItem("user_id");
      if (!userId) {
        this.$message({
          showClose: true,
          message: "登录已过期，请重新登录",
          type: "error",
        });
        this.$router.push("/login");
        return;
      }

      //创建formdata对象
      const formData = new FormData();
      // 添加用户ID到表单数据
      formData.append("user_id", userId);
      formData.append("name", this.form.name);
      formData.append("species", this.form.species);
      formData.append("desc", this.form.desc);
      formData.append("save_path", this.form.save_path);
      formData.append("figure_path", this.form.figure_path);

      this.$refs[formName].validate((valid) => {
        if (valid) {
          this.$message("正在创建中，请稍等……");
          this.dialogFormVisible = false;

          // 使用 api.js 中的 post 方法
          post("create_item", formData)
            .then((response) => {
              this.$message({
                showClose: true,
                message: "项目已创建",
                type: "success",
              });
              console.log(response.data);
              localStorage.setItem('item_id', response.data.item_id);

              this.$router.push("/sdcs");
            })
            .catch((error) => {
              this.$message({
                showClose: true,
                message: "创建失败：" + error.message,
                type: "error",
              });
            });
        } else {
          console.log("error submit!!");
          return false;
        }
      });
    },

    // 添加跳转方法
    jumpToTutorial(section) {
      // 修改路径为根路径
      window.open(`/tutorial.html#${section}`, "_blank");
    },

    async showItems() {
      try {
        const response = await get("get_user_items");
        if (response.data.status === "success") {
          this.items = response.data.items;
        } else {
          this.$message.error(response.data.message || "获取项目列表失败");
        }
      } catch (error) {
        this.$message.error("获取项目列表失败：" + error.message);
      }
    },

    getSpeciesTagType(species) {
      const types = {
        human: "success",
        mouse: "warning",
        other_species: "info",
      };
      return types[species] || "info";
    },

    async viewProject(item) {
      localStorage.setItem('item_id', item.id);
      this.projectName = item.name;
      this.taskDialogVisible = true;
      this.taskLoading = true;
      try {
        const res = await get('get_item_tasks');
        if (res.data.status === "success") {
          // 更新树状图数据
          this.treeData = res.data.tree;
          // 更新表格数据
          this.taskList = res.data.tasks;
          this.$nextTick(() => {
            this.renderTaskTree();
          });
        } else {
          this.$message.error(res.data.message || "获取任务失败");
        }
      } catch (e) {
        this.$message.error("获取任务失败：" + e.message);
      }
      this.taskLoading = false;
    },

    renderTaskTree() {
      if (!this.treeData) return;
      const treeData = [this.treeData];
      const chartDom = this.$refs.treeChart;
      if (!chartDom) return;
      const myChart = echarts.init(chartDom);

      const option = {
        backgroundColor: '#f8fafc',
        tooltip: {
          trigger: 'item',
          triggerOn: 'mousemove',
          backgroundColor: '#fff',
          borderColor: '#409EFF',
          borderWidth: 1,
          textStyle: {
            color: '#333',
            fontSize: 14
          },
          padding: 10,
          extraCssText: 'box-shadow: 0 2px 8px rgba(64,158,255,0.15);'
        },
        series: [
          {
            type: 'tree',
            data: treeData,
            top: '5%',
            left: '10%',
            bottom: '5%',
            right: '25%',
            symbol: 'circle',
            symbolSize: 18,
            edgeShape: 'polyline',
            edgeForkPosition: '63%',
            roam: true,
            initialTreeDepth: 3,
            nodePadding: 10, // 节点之间的垂直间距，默认7，调大可拉开条目
            lineStyle: {
              color: '#409EFF',
              width: 2,
              curveness: 0.25,
              shadowColor: 'rgba(64,158,255,0.2)',
              shadowBlur: 6
            },
            itemStyle: {
              color: '#fff',
              borderColor: '#409EFF',
              borderWidth: 2,
              shadowColor: 'rgba(64,158,255,0.15)',
              shadowBlur: 8
            },
            label: {
              backgroundColor: '#fff',
              borderColor: '#409EFF',
              borderWidth: 1,
              borderRadius: 6,
              padding: [6, 12],
              color: '#333',
              fontWeight: 'bold',
              fontSize: 15,
              position: 'left',
              verticalAlign: 'middle',
              align: 'right',
              shadowColor: 'rgba(64,158,255,0.08)',
              shadowBlur: 2,
            },
            leaves: {
              label: {
                backgroundColor: '#e6f7ff',
                borderColor: '#91d5ff',
                borderWidth: 1,
                borderRadius: 6,
                padding: [2, 10],
                color: '#409EFF',
                fontWeight: 'normal',
                fontSize: 12,
                position: 'right',
                verticalAlign: 'middle',
                align: 'left',
                minMargin:22,
                height: 22,
                lineHeight: 22 // 叶子节点的行高，进一步拉开条目
              }
            },
            expandAndCollapse: true,
            animationDuration: 600,
            animationDurationUpdate: 800
          }
        ]
      };
      myChart.setOption(option);
      window.addEventListener('resize', () => {
        myChart.resize();
      });
    },

    // 预留删除方法
    deleteProject(item) {
      this.$confirm('确认删除该项目吗？此操作不可恢复', '提示', {
        confirmButtonText: '确定',
        cancelButtonText: '取消',
        type: 'warning'
      }).then(async () => {
        localStorage.setItem('item_id', item.id);
        try {
          const response = await post('delete_item');
          if (response.data.status === 'success') {
            this.$message({
              type: 'success',
              message: '删除成功!'
            });
            // 重新加载项目列表
            this.showItems();
          } else {
            this.$message.error(response.data.message || '删除失败');
          }
        } catch (error) {
          this.$message.error('删除失败：' + error.message);
        }
      }).catch(() => {
        this.$message({
          type: 'info',
          message: '已取消删除'
        });
      });
    },

    handleTaskDelete(row) {
      this.$confirm('确认删除该任务记录吗？此操作不可恢复', '提示', {
        confirmButtonText: '确定',
        cancelButtonText: '取消',
        type: 'warning'
      }).then(async () => {
        try {
          const response = await post('delete_task', {
            task_id: row.task_id
          });
          if (response.data.status === 'success') {
            this.$message({
              type: 'success',
              message: '删除成功!'
            });
            // 重新加载任务列表和树状图
            const res = await get('get_item_tasks');
            if (res.data.status === "success") {
              this.taskList = res.data.tasks;
              this.treeData = res.data.tree;  // 更新树状图数据
              this.$nextTick(() => {
                this.renderTaskTree();  // 重新渲染树状图
              });
            }
          } else {
            this.$message.error(response.data.message || '删除失败');
          }
        } catch (error) {
          this.$message.error('删除失败：' + error.message);
        }
      }).catch(() => {
        this.$message({
          type: 'info',
          message: '已取消删除'
        });
      });
    },

    handleCardClick(item) {
      localStorage.setItem('item_id', item.id);
      this.$router.push('/jsjg');  // 假设可视化结果页面的路由是 /jsjg
    },

    handleTaskView(row) {
     
      if (row.type) {
        // 通过路由导航到jsjg页面，并传递目标tab参数
        this.$router.push({
          path: '/jsjg',
          query: { 
            type: row.type,
            taskId: row.task_id           // 可以传递任务ID，以便在目标页面加载对应数据
          }
        })
      } else {
        this.$message.warning('未知的任务类型')
      }
    },
  },
};
</script>

<style scoped>
.block {
  margin-top: 0px;
  display: flex;
  justify-content: space-between;
}

.insert {
  width: 10%;
}

.container {
  width: 80%;
  align-items: center;
  display: flex;
  justify-content: space-between;
  border: 4px solid deepskyblue;
  border-radius: 10px;
  /* 添加圆角 */
}

.table_head {
  height: 15%;
  width: 100%;
  display: flex;
  justify-content: center;
  align-items: center;
  color: #f4fcfe;
  font-family: SanJiKaiShu;
  font-size: 23px;
  /* background-color: #007bff; */
  background-color: deepskyblue;
}

.panel {
  width: 100%;
  height: 100%;
}

.rounded-card {
  border-radius: 10px;
  /* 添加圆角 */
}

.dialog-footer {
  display: flex;
  justify-content: flex-end;
}

.el-row {
  margin-bottom: 20px;

  &:last-child {
    margin-bottom: 0;
  }
}

.el-col {
  border-radius: 4px;
}

.bg-purple-dark {
  background: #99a9bf;
}

.bg-purple {
  background: #d3dce6;
}

.grid-content {
  border-radius: 4px;
  min-height: 36px;
}

.row-bg {
  padding: 10px 0;
  background-color: #f9fafc;
}

.projects-header {
  display: flex;
  align-items: center;
  margin-bottom: 20px;
  gap: 30px;
}

.button-group {
  display: flex;
  align-items: center;
  gap: 0px;
}

.help-button {
  padding: 2px;
  min-width: 24px;
  height: 24px;
}

.projects-grid {
  display: grid;
  grid-template-columns: repeat(auto-fill, minmax(280px, 1fr));
  gap: 20px;
  margin-top: 20px;
}

.project-card {
  transition: all 0.3s;
  cursor: pointer;  /* 添加鼠标指针样式 */
}

.project-card:hover {
  transform: translateY(-5px);
  box-shadow: 0 2px 12px 0 rgba(0,0,0,.1);  /* 添加悬停阴影效果 */
}

.project-info {
  display: flex;
  flex-direction: column;
  gap: 4px;
}

.project-time {
  font-size: 12px;
  color: #999;
}

.project-header {
  display: flex;
  justify-content: space-between;
  align-items: flex-start;
}

.project-name {
  font-weight: bold;
  font-size: 16px;
}

.project-content {
  min-height: 100px;
  display: flex;
  flex-direction: column;
  justify-content: space-between;
}

.project-desc {
  color: #666;
  margin: 10px 0;
  display: -webkit-box;
  -webkit-line-clamp: 2;
  -webkit-box-orient: vertical;
  overflow: hidden;
}

.project-actions {
  display: flex;
  justify-content: flex-end;
  margin-top: 10px;
  gap: 0px;  /* 添加按钮间距 */
}

.action-button {
  padding: 0px 10;  /* 减小内边距 */
  min-width: unset;  /* 移除最小宽度限制 */
}
.action-button,
.el-button--primary,
.el-button--danger {
  border-radius: 20px !important;
  min-width: 30px;
}
</style>
