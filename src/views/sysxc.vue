<template>
  <div class="lab-scene">
    <el-card>
      <div class="select-container">
        <h2>养殖仓编号</h2>
        <el-select v-model="value" class="m-2 select" placeholder="Select">
          <el-option
            v-for="item in options"
            :key="item.value"
            :label="item.label"
            :value="item.value"
          />
        </el-select>
      </div>

      <el-divider></el-divider>
      <div class="camera-control">
        <h3>摄像头操作</h3>
        <div class="xunjian">
          <el-button
            type="primary"
            style="width: 100%; margin-bottom: 10px"
            @click="this.operate('auto')"
            >自动巡检</el-button
          >
          <el-button
            type="primary"
            style="width: 100%; margin-bottom: 10px"
            @click="this.operate('continue')"
            >继续巡检</el-button
          >
        </div>

        <div class="directions">
          <el-button type="primary" @click="operate('move up')"
            >上</el-button
          >
          <el-button type="primary" @click="operate('move down')"
            >下</el-button
          >
          <el-button type="primary" @click="operate('move left')"
            >左</el-button
          >
          <el-button type="primary" @click="operate('move right')"
            >右</el-button
          >
          <el-button type="primary" @click="operate('pause')"
            >暂停</el-button
          >
        </div>
        <div class="input">
          <span style="font-weight: 700;margin-top:5px;">移动至</span>
          <el-input
            v-model="inputValue"
            style="width: 100px"
            placeholder="请输入内容"
          >
          </el-input
          ><el-button @click="operate(`book ${inputValue}`)">发送命令</el-button>
        </div>
      </div>
      <el-divider></el-divider>
      <h3>其它操作</h3>
      <el-button type="success" @click="operate('')">自动投食</el-button>
      <el-button type="success" @click="operate('')">自动喂水</el-button>
    </el-card>
    <el-card class="right-panel">
      <h2>实时显示摄像头画面</h2>
      <el-divider></el-divider>
      摄像头当前位置：{{ pos }}
      <div class="camera-frame">
        <!-- 摄像头画面展示区域 -->
        <img
          v-if="value"
          style="width: 100%"
          src="api/video_feed"
          ref="video"
        />
      </div>
    </el-card>
  </div>
</template>

  <script>
import { send_move } from "../request/api.js";
import { ElMessage } from "element-plus";
import { ElInput, ElButton } from "element-plus";
import axios from 'axios'

export default {
  components: {},
  data() {
    return {
      inputValue: "",
      value: "",
      options: [
        {
          value: "1",
          label: "1",
        },
        {
          value: "2",
          label: "2",
        },
        {
          value: "3",
          label: "3",
        },
        {
          value: "4",
          label: "4",
        },
        {
          value: "5",
          label: "5",
        },
      ],
      pos:null,
      timer:null
    };
  },
  methods: {
    operate(message) {
      this.inputValue = "";
      send_move(message)
        .then((res) => {
          ElMessage({
            message: `发送成功:${JSON.stringify(res)}`,
            type: "success",
          });
        })
        .catch((err) => {
          ElMessage({
            message: "发送失败" + JSON.stringify(err),
            type: "error",
          });
        });
    },
    getPos() {
      axios.get('/api/check_pos')
        .then(response => {
          if (response.data && response.data.pos) {
            console.log('response.data',response.data)
            this.pos = response.data.pos
          }
        })
        .catch(error => {
          console.error(error)
        })
    }
  },
  created(){
    this.getPos() // initial fetch
    // this.timer = setInterval(this.getPos, 3000) // fetch every 5 seconds
  },
  // beforeDestroy() {
  //   clearInterval(this.timer) // clear the timer when the component is unmounted
  // },
  watch: {
    value() {
      ElMessage({
        message: "切换成功",
        type: "success",
      });
    },
  },
};
</script>
  
  
  <style>
.lab-scene {
  margin-top: 0px;
  display: flex;
  justify-content: space-between;
  /* border: 1px solid #ccc; */
  padding: 10px;
}

.left-panel {
  width: 30%;
  padding: 10px;
  border-right: 1px solid #ccc;
}

.camera-control {
  display: flex;
  flex-direction: column;
  margin-bottom: 10px;
}

.right-panel {
  width: 70%;
  padding: 10px;
}

.camera-frame {
  border: 1px solid #ccc;
  height: 350px;
}
</style>

  <style>
.el-main {
  --el-main-padding: 0px 20px !important;
}
.select-container {
  display: flex;
  align-items: center;
}

.select {
  margin-left: 10px;
  width: 100px;
}

.directions {
  display: flex;
  justify-content: space-between;
}
.input {
  margin-top: 10px;
  display: flex;
  justify-content: space-between;
}
.xunjian {
  display: flex;
  flex-direction: row;
}
.el-button+.el-button {
    margin-left: 5px;
}
</style>
  