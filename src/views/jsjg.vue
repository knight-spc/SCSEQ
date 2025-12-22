<template>
  <el-tabs v-model="activeName" @tab-click="handleClickTab">
    <el-tab-pane label="单细胞测序结果" name="first">
      <n-split direction="horizontal" style="height: 100%" default-size="500px">
        <template #1>
          Dimplot<el-tooltip
            class="box-item"
            effect="dark"
            content="介绍手册"
            placement="top"
          >
            <el-button
              type=""
              text
              size="small"
              @click="jumpToTutorial('annotation-view')"
            >
              <el-icon>
                <QuestionFilled />
              </el-icon>
            </el-button>
          </el-tooltip>
          <div>
            <el-button
              type="primary"
              :loading="cluster_figure.if_load"
              @click="getUmapData"
              >点击查看umap图</el-button
            >
            <el-tree-select
              v-model="annotation.value"
              :data="annotation.models"
              :render-after-expand="false"
              :placeholder="annotation.placeholder"
              @focus="fetchAnnotationModels"
              @change="handleSelectChange_annotation"
              style="width: 220px"
            />
            <el-dialog
              v-model="annotation.dialogVisible"
              title="上传数据"
              width="60%"
              height="100%"
            >
              <el-row :gutter="30">
                <el-col :span="12">
                  <el-form :model="annotation" ref="annotation">
                    <el-form-item
                      label="模型名称"
                      prop="modelName"
                      :rules="[
                        { required: true, message: '模型名称不能为空' },
                        {
                          pattern: /^[a-zA-Z0-9_]+$/,
                          message: '模型名称只能包含字母、数字和下划线',
                        },
                      ]"
                    >
                      <el-input
                        v-model="annotation.modelName"
                        autocomplete="off"
                      ></el-input>
                    </el-form-item>
                  </el-form>
                  <div class="block">
                    <div class="insert"></div>
                    <div class="container">
                      <div class="panel">
                        <div class="table_head">请在这里上传参考数据集</div>
                        <el-upload
                          :headers="uploadHeaders"
                          class="upload-demo"
                          drag
                          name="zip"
                          action="/api/upload_reference"
                          accept=".zip"
                          :on-success="handleReferenceDataSuccess"
                          :on-error="() => $message.error('上传')"
                          style="
                            flex-direction: column;
                            display: flex;
                            justify-content: center;
                          "
                        >
                          <el-icon class="el-icon--upload"><upload-filled /></el-icon>
                          <div class="el-upload__text">
                            Drop file here or <em>click to upload</em>
                          </div>
                          <template #tip>
                            <div class="el-upload__tip" style="text-align: center">
                              注意：上传前请保证文件格式正确
                            </div>
                          </template>
                        </el-upload>
                      </div>
                    </div>
                    <div class="insert"></div>
                  </div>
                </el-col>
                <el-col :span="12">
                  <p>上传文件的格式要求：<br />请将您的注释放入Annotation列中</p>

                  <el-row>
                    <div class="grid-content">
                      <el-image 
                        :src="reference_example" 
                        fit="contain"
                        style="max-width: 100%; max-height: 200px; display: block; margin: 0 auto;"
                      ></el-image>
                    </div>
                  </el-row>
                </el-col>
              </el-row>

              <div slot="footer" class="dialog-footer">
                <el-button @click="annotation.dialogVisible = false">取 消</el-button>
                <el-button type="primary" @click="submitReferenceData">确 定</el-button>
              </div>
            </el-dialog>
            <!-- 用于显示ECharts图表的div，设置了宽度和高度 -->
            <div ref="chart" style="width: 500px; height: 500px"></div>
          </div>
          <el-button type="primary" :loading="marker.if_load" @click="CalculateMarkerGene"
            >计算Marker基因</el-button
          >
          <el-button type="primary" :loading="marker.if_load" @click="LargeModelAnnotation"
            >大模型辅助注释</el-button
          >

    <!-- 大模型结果弹窗 -->
    <el-dialog v-model="annotation.LargeModelDialogVisible" title="大模型辅助预测结果" width="80%" :modal="true" :close-on-click-modal="false" style="max-width:1200px;">
      
      <el-table :data="annotation.tableData" style="width: 100%; margin-top: 20px;">
        <el-table-column prop="origin" label="原注释结果" width="200"/>
        <el-table-column prop="large" label="大模型注释结果" width="200"/>
        <el-table-column prop="reason" label="推断依据"  />
      </el-table>
    </el-dialog>

        </template>
        <template #2>
          <n-split direction="horizontal" style="height: 100%">
            <template #1>
              QC小提琴图<el-tooltip
                class="box-item"
                effect="dark"
                content="介绍手册"
                placement="top"
              >
                <el-button
                  type=""
                  text
                  size="small"
                  @click="jumpToTutorial('quality-control')"
                >
                  <el-icon>
                    <QuestionFilled />
                  </el-icon>
                </el-button>
              </el-tooltip>
              <!-- <div style="width=10px,height: 10px">
                    <el-image :src="url" :fit="contain" style="width=100px,height: 100px"></el-image>
                  </div> -->
              <div class="demo-image__lazy">
                <!-- <div class="insert"></div> -->
                <el-scrollbar ref="el_scrollbar" style="height: 600px; overflow-y: auto">
                  <el-image
                    v-for="path in imgPaths"
                    :key="path"
                    :src="require('../assets/result/' + path)"
                    style="height: 350px; margin: 10px"
                    fit="contain"
                    lazy
                    :scroll-container="scrollContainer"
                    :preview-src-list="[require('../assets/result/' + path)]"
                  >
                    <template #toolbar="{ actions, reset }">
                      <el-icon @click="actions('zoomOut')">
                        <ZoomOut />
                      </el-icon>
                      <el-icon
                        @click="
                          actions('zoomIn', { enableTransition: false, zoomRate: 2 })
                        "
                      >
                        <ZoomIn />
                      </el-icon>
                      <el-icon
                        @click="
                          actions('clockwise', {
                            rotateDeg: 180,
                            enableTransition: false,
                          })
                        "
                      >
                        <RefreshRight />
                      </el-icon>
                      <el-icon @click="actions('anticlockwise')">
                        <RefreshLeft />
                      </el-icon>
                      <el-icon @click="reset">
                        <Refresh />
                      </el-icon>
                      <el-icon
                        @click="downloadImage(require('../assets/result/' + path))"
                      >
                        <Download />
                      </el-icon>
                    </template>
                    <template #error>
                      <div class="image-slot">
                        <el-icon><icon-picture /></el-icon>
                      </div>
                    </template>
                  </el-image>
                </el-scrollbar>
                <!-- <div class="insert"></div> -->
              </div>
            </template>
            <template #2>
              基于自动注释的细胞比例图<el-tooltip
                class="box-item"
                effect="dark"
                content="介绍手册"
                placement="top"
              >
                <el-button
                  type=""
                  text
                  size="small"
                  @click="jumpToTutorial('proportion-view')"
                >
                  <el-icon>
                    <QuestionFilled />
                  </el-icon>
                </el-button>
              </el-tooltip>

              <el-scrollbar ref="el_scrollbar" style="height: 600px; overflow-y: auto">
                <div class="znjy-container">
                  <div class="controls">
                    <el-select
                      v-model="scale_diagram.selectedView"
                      placeholder="选择视图"
                      @change="handleViewChange"
                      style="width: 150px"
                    >
                      <el-option
                        v-for="item in scale_diagram.viewOptions"
                        :key="item.value"
                        :label="item.label"
                        :value="item.value"
                      >
                      </el-option>
                    </el-select>


                    <el-select-v2
                      v-model="scale_diagram.categoryMode"
                      filterable
                      :options="scale_diagram.models"
                      @change="applyCategoryMode"
                      placeholder="Please select"
                      style="width: 200px"
                    />
                  </div>
                  <div class="chart-container">
                    <div
                      ref="chart_scale_diagram"
                      id="chart_scale_diagram"
                      style="width: 400px; height: 400px"
                    ></div>
                  </div>
                </div>
              </el-scrollbar>
            </template>
          </n-split>
        </template>
      </n-split>
    </el-tab-pane>

    <el-tab-pane label="基因富集分析" name="second">
      <n-split direction="horizontal" style="height: 100%" default-size="470px">
        <template #1>
          <n-split direction="vertical" style="height: 100%; width: 100%">
            <template #1>
              Table
              <el-tooltip
                class="box-item"
                effect="dark"
                content="介绍手册"
                placement="top"
              >
                <el-button
                  type=""
                  text
                  size="small"
                  @click="jumpToTutorial('table-analysis')"
                >
                  <el-icon>
                    <QuestionFilled />
                  </el-icon>
                </el-button>
              </el-tooltip>
              <el-table
                :data="marker.tableData"
                max-height="250"
                lazy
                fit="false"
                stripe="true"
                border="true"
                resizable
                show-overflow-tooltip="true"
              >
                <el-table-column fixed prop="gene" label="gene" width="80">
                </el-table-column>
                <el-table-column prop="p_val" label="p_val" width="70"> </el-table-column>
                <el-table-column prop="avg_log2FC" label="avg_log2FC" width="80">
                </el-table-column>
                <el-table-column prop="pct.1" label="pct.1" width="60"> </el-table-column>
                <el-table-column prop="pct.2" label="pct.2" width="60"> </el-table-column>
                <el-table-column prop="p_val_adj" label="p_val_adj" width="60">
                </el-table-column>
                <el-table-column prop="cluster" label="cluster" width="60">
                </el-table-column>
              </el-table>
              <el-row>
                <el-select
                  v-model="marker.num"
                  filterable
                  allow-create
                  default-first-option
                  style="width: 30%"
                  @change="handleSelectChange_marker"
                  placeholder="请选择每组标记基因的数目"
                >
                  <el-option
                    v-for="item in marker.num_list"
                    :key="item.value"
                    :label="item.label"
                    :value="item.value"
                  >
                  </el-option>
                </el-select>
                <el-button type="primary" class="mt-4" @click="exportToCSV"
                  >Download as CSV</el-button
                >
              </el-row>
            </template>
            <template #2>
              DotPlot
              <el-tooltip
                class="box-item"
                effect="dark"
                content="介绍手册"
                placement="top"
              >
                <el-button
                  type=""
                  text
                  size="small"
                  @click="jumpToTutorial('dotplot-analysis')"
                >
                  <el-icon>
                    <QuestionFilled />
                  </el-icon>
                </el-button>
              </el-tooltip>
              <el-select
                v-model="marker.dotplot_num"
                filterable
                allow-create
                default-first-option
                style="width: 30%"
                @change="getDotPlotFigure"
                placeholder="请选择每组标记基因的数目"
              >
                <el-option
                  v-for="item in marker.num_list"
                  :key="item.value"
                  :label="item.label"
                  :value="item.value"
                >
                </el-option>
              </el-select>
              <div style="position: relative; width: 1150px">
                <!-- X轴图层 -->
                <div
                  ref="xAxisChart"
                  style="
                    position: relative;
                    width: 1000px;
                    height: 80px;
                    margin-left: 150px;
                  "
                ></div>

                <!-- 可滚动区域 -->
                <el-scrollbar ref="el_scrollbar" style="height: 500px; margin-top: 5px">
                  <div style="position: relative; height: 1000px; width: 1150px">
                    <!-- 合并后的图表 -->
                    <div
                      ref="chart_dotplot"
                      style="
                        position: absolute;
                        top: 0;
                        left: 0;
                        width: 1150px;
                        height: 1000px;
                      "
                    ></div>
                  </div>
                </el-scrollbar>
              </div>
              <!-- <div ref="chart_dotplot" style="width: 100%; height: 300px;"></div> -->
            </template>
          </n-split>
        </template>
        <template #2>
          <n-split direction="horizontal" style="height: 100%; width: 100%">
            <template #1>
              <!-- <n-split direction="vertical" style="height: 100%; width: 100%">
                <template #1> -->
              Go富集
              <el-tooltip
                class="box-item"
                effect="dark"
                content="介绍手册"
                placement="top"
              >
                <el-button
                  type=""
                  text
                  size="small"
                  @click="jumpToTutorial('go-enrichment')"
                >
                  <el-icon>
                    <QuestionFilled />
                  </el-icon>
                </el-button>
              </el-tooltip>
              <el-select
                v-model="go.choose_cluster"
                filterable
                lazy
                placeholder="Select"
                @change="handleSelectChange_go"
                style="width: 350px"
              >
                <el-option
                  v-for="item in go.clusters"
                  :key="item.value"
                  :label="item.label"
                  :value="item.value"
                />
              </el-select>
              <div ref="chart_bar" style="width: 500px; height: 500px"></div>
              <!-- </template> -->
              <!-- <template #2>
                  panel4
                  <div><el-image :src="url" :fit="contain"></el-image></div>
                </template> -->
              <!-- </n-split> -->
            </template>
            <template #2>
              <n-split direction="vertical" style="height: 100%">
                <template #1>
                  QC小提琴图
                  <el-tooltip
                    class="box-item"
                    effect="dark"
                    content="介绍手册"
                    placement="top"
                  >
                    <el-button
                      type=""
                      text
                      size="small"
                      @click="jumpToTutorial('qc-violin')"
                    >
                      <el-icon>
                        <QuestionFilled />
                      </el-icon>
                    </el-button>
                  </el-tooltip>

                  <el-button
                    type="primary"
                    style="margin-left: 16px"
                    @click="SignificanceQC.drawer = true"
                  >
                    查看显著性结果
                  </el-button>

                  <el-drawer
                    v-model="SignificanceQC.drawer"
                    title="I am the title"
                    :with-header="false"
                  >
                    <span>选择两个细胞类型</span>
                    <el-checkbox-group
                      v-model="SignificanceQC.checkedCells"
                      :min="0"
                      :max="2"
                    >
                      <el-checkbox
                        v-for="cell in SignificanceQC.citiescellOptions"
                        :label="cell"
                        :key="cell"
                        >{{ cell }}</el-checkbox
                      >
                    </el-checkbox-group>
                    <el-button
                      type="primary"
                      @click="submitSelection_Significance"
                      :loading="SignificanceQC.loading"
                      >提交选择</el-button
                    >
                    <div v-if="SignificanceQC.imagePath" class="image-container">
                      <el-image
                        v-for="path in SignificanceQC.imagePath"
                        :key="path"
                        :src="require('../assets/result/' + path)"
                        style="height: 350px; margin: 10px"
                        fit="contain"
                        lazy
                        :preview-src-list="[require('../assets/result/' + path)]"
                      >
                        <template #toolbar="{ actions, reset }">
                          <el-icon @click="actions('zoomOut')">
                            <ZoomOut />
                          </el-icon>
                          <el-icon
                            @click="
                              actions('zoomIn', { enableTransition: false, zoomRate: 2 })
                            "
                          >
                            <ZoomIn />
                          </el-icon>
                          <el-icon
                            @click="
                              actions('clockwise', {
                                rotateDeg: 180,
                                enableTransition: false,
                              })
                            "
                          >
                            <RefreshRight />
                          </el-icon>
                          <el-icon @click="actions('anticlockwise')">
                            <RefreshLeft />
                          </el-icon>
                          <el-icon @click="reset">
                            <Refresh />
                          </el-icon>
                          <el-icon
                            @click="downloadImage(require('../assets/result/' + path))"
                          >
                            <Download />
                          </el-icon>
                        </template>
                        <template #error>
                          <div class="image-slot">
                            <el-icon><icon-picture /></el-icon>
                          </div>
                        </template>
                      </el-image>
                    </div>
                  </el-drawer>

                  <!-- <div style="width=10px,height: 10px">
                    <el-image :src="url" :fit="contain" style="width=100px,height: 100px"></el-image>
                  </div> -->
                  <div class="demo-image__lazy">
                    <!-- <div class="insert"></div> -->
                    <el-scrollbar
                      ref="el_scrollbar"
                      style="height: 300px; overflow-y: auto"
                    >
                      <el-image
                        v-for="path in marker.imgPaths"
                        :key="path"
                        :src="require('../assets/result/' + path)"
                        style="height: 350px; margin: 10px"
                        fit="contain"
                        lazy
                        :scroll-container="scrollContainer"
                        :preview-src-list="[require('../assets/result/' + path)]"
                      >
                        <template #toolbar="{ actions, reset }">
                          <el-icon @click="actions('zoomOut')">
                            <ZoomOut />
                          </el-icon>
                          <el-icon
                            @click="
                              actions('zoomIn', { enableTransition: false, zoomRate: 2 })
                            "
                          >
                            <ZoomIn />
                          </el-icon>
                          <el-icon
                            @click="
                              actions('clockwise', {
                                rotateDeg: 180,
                                enableTransition: false,
                              })
                            "
                          >
                            <RefreshRight />
                          </el-icon>
                          <el-icon @click="actions('anticlockwise')">
                            <RefreshLeft />
                          </el-icon>
                          <el-icon @click="reset">
                            <Refresh />
                          </el-icon>
                          <el-icon
                            @click="downloadImage(require('../assets/result/' + path))"
                          >
                            <Download />
                          </el-icon>
                        </template>
                        <template #error>
                          <div class="image-slot">
                            <el-icon><icon-picture /></el-icon>
                          </div>
                        </template>
                      </el-image>
                    </el-scrollbar>
                    <!-- <div class="insert"></div> -->
                  </div>
                </template>

              </n-split>
            </template>
          </n-split>
        </template>
      </n-split>
    </el-tab-pane>

  </el-tabs>
</template>
<script>
import axios from "axios";
import { ElMessage, ElMessageBox } from "element-plus";
import { Eleme } from "@element-plus/icons-vue";
import { UploadFilled } from "@element-plus/icons-vue";
import * as echarts from "echarts"; // 引入ECharts库
import Papa from "papaparse"; // 引入PapaParse库，用于解析CSV文件
import { SomeComponent } from "naive-ui"; // 替换为实际使用的组件
import { ElNotification } from "element-plus";

import { NInputNumber, NSplit } from "naive-ui";
import { ref } from "vue";
const split = ref(0.8);
import { markRaw } from "vue";
import { get, post } from "@/request/api";

const chartInstance = ref(null);

import {
  Download,
  Refresh,
  RefreshLeft,
  RefreshRight,
  ZoomIn,
  ZoomOut,
  Picture as IconPicture,
} from "@element-plus/icons-vue";

export default {
  components: {
    UploadFilled,
    IconPicture,
    Download,
    Refresh,
    RefreshLeft,
    RefreshRight,
    ZoomIn,
    ZoomOut,
  },
  data() {
    return {
      activeName: "first",
      url: require("../assets/images/scseq2.png"),
      imgPaths: [],
      if_load: false,
      split: 0,
      task_id:'',
      cluster_figure: {
        if_load: false,

        chart: null, // ECharts实例，用于后续操作图表
        clustered_data: [], // 存储聚类后的数据
        // clustered_label: [], // 存储聚类后的数据标签
        clustered_distribution: [], // 存储聚类数据描述
        clusters: [], // 存储聚类结果，每个聚类包含标签、中心点和数据点
        cluster_colors: [], // 为不同聚类分配的颜色
        legend_label: [], // 图例数据
      },
      annotation: {
        value: "",
        models: [
        {vaule:"大模型",
        label:"大模型",
        children: [
              {
                value: "AICelltype",
                label: "AICelltype",
              }]
        },
          {
            value: "SingleR",
            label: "SingleR",
            children: [
              {
                value: "HumanPrimaryCellAtlasData",
                label: "HumanPrimaryCellAtlasData",
              },
              {
                value: "BlueprintEncodeData",
                label: "BlueprintEncodeData",
              },
              {
                value: "DatabaseImmuneCellExpressionData",
                label: "DatabaseImmuneCellExpressionData",
              },
              {
                value: "NovershternHematopoieticData",
                label: "NovershternHematopoieticData",
              },
              {
                value: "MonacoImmuneData",
                label: "MonacoImmuneData",
              },
              {
                value: "MouseRNAseqData",
                label: "MouseRNAseqData",
              },
              {
                value: "ImmGenData",
                label: "ImmGenData",
              },
            ],
          },
        ],
        dialogVisible: false,
        reference_path: "",
        currunt_model: "seurat_clusters",
        placeholder: "请选择参考数据集",
        modelName: "",
        save_path: "",
        tableData: [], // 存储表中的数据
        LargeModelDialogVisible:false,
      },
      scale_diagram: {
        myChart: null,
        selectedView: "percentage",
        viewOptions: [
          { value: "percentage", label: "百分比视图" },
          { value: "absolute", label: "绝对值视图" },
        ],
        // 这些数据将从后端获取
        rawData: [],
        totalData: [],
        categories: [],
        models: [],
        categoryMode: "orig.ident",
        seriesNames: "",
        serie:"", // 用于绘制比例图时确定cluster还是注释
      },

      marker: {
        chart: null, // ECharts实例，用于后续操作图表
        tableData: [], // 存储表中的数据
        num_list: [
          {
            value: "1",
            label: "1",
          },
          {
            value: "2",
            label: "2",
          },
          {
            value: "5",
            label: "5",
          },
          {
            value: "10",
            label: "10",
          },
          {
            value: "20",
            label: "20",
          },
        ], // 基因数目下拉框选项
        num: "",
        gene_list: [], // 特征基因下拉框选项
        gene: "",
        clustered_data: [], // 存储聚类后的数据
        clusters: [], // 存储聚类结果，每个聚类包含标签、中心点和数据点
        gene_expression: [], // 基因表达数据
        imgPaths: [], // 存储Vln图片路径
        dotplot_data: [], // 存储dotplot数据
        dotplot_num: "", // 存储dotplot数据
        if_load: false,
      },
      SignificanceQC: {
        checkedCells: [],
        citiescellOptions: [],
        loading: false,
        imagePath: [],
        drawer: false, // 用于显示抽屉
      },

      go: {
        chart: null,
        value: "",
        choose_cluster: "",
        bar_data: [],
        clusters: [], // 存储下拉框选项，用于存储细胞群序号
      },
      reference_example: require("../assets/images/reference_example.png"),
    };
  },
  created() {
      // 从路由参数中获取目标tab
      const type = this.$route.query.type
      // 任务类型与tab名称的映射
      const taskTypeToTab = {
        'Annotation': 'first',          // QC任务对应单细胞测序结果tab
        'Marker': 'second',          // Marker任务对应基因富集分析tab
        'Feature': 'second',          // Marker任务对应基因富集分析tab
        'Go': 'second',          // Marker任务对应基因富集分析tab
        'CellChat': 'third',         // CellChat任务对应细胞通讯分析tab
        'InferCNV': 'fourth',        // InferCNV任务对应infercnv tab
        'Pseudotime': 'fifth',        // Pseudotime任务对应拟时序分析tab
        'Pseudotime_gene': 'fifth',        // Pseudotime任务对应拟时序分析tab
        'TCGA': 'sixth',        // Pseudotime任务对应拟时序分析tab
        'TCGA_genes': 'seventh',        // Pseudotime任务对应拟时序分析tab
      }
      const activeTab = taskTypeToTab[type]

      if (activeTab) {
        this.activeName = activeTab
        // 可以根据taskId加载对应的数据
        this.task_id = this.$route.query.taskId

        if (this.task_id) {
          this.loadTaskData(type,this.task_id)
        }
      }
    },

  computed: {
    uploadHeaders() {
      return {
        "User-Id": localStorage?.getItem("user_id") || "",
        "Item-Id": localStorage?.getItem("item_id") || "",
      };
    },
  },
  methods: {
    async loadTaskData(type, taskId) {
      switch (type) {
        case 'Annotation':
          ElMessage.success("即将展示Annotation数据");
          this.getUmapData();
          break;
        case 'Marker':
          ElMessage.success("即将展示Marker数据");
          this.handleSelectChange_marker(5);
          this.getDotPlotFigure(5)
          break;
        case 'Go':
          this.fetchGoFigure();
          break;
      }

    },

    handleClickTab(tab, event) {
      console.log(tab, event);
      if (tab.props.name === "third") {
        // 如果切换到细胞通讯分析标签页
        this.$nextTick(() => {
          // 等待DOM更新后再获取引用
          this.cellchat.scrollContainer = this.$refs.cellchat_scrollbar.$el.wrap;
          // 如果热点图路径为空，则获取热点图数据
          if (this.cellchat.heatmapPaths.length === 0) {
            this.getCellChatHeatmaps();
            this.initCircleChart();
          }
        });
      }
      // 当切换到 infercnv 标签页时
      if (tab.props.name === "fourth" || tab.props.name === "fifth") {
        // this.infercnv.loading = true; // 显示加载状态
        this.getCellList()
          .then(() => {
            this.infercnv.loading = false; // 数据加载完成后关闭加载状态
          })
          .catch((error) => {
            console.error("获取细胞列表失败:", error);
            this.$message.error("获取细胞列表失败");
            this.infercnv.loading = false;
          });
      }
    },

    beforeDestroy() {
      // 在组件销毁前，销毁ECharts实例以释放资源
      if (this.chart) {
        this.chart.dispose();
      }
    },
    // 获取UMAP图数据
    getUmapData() {
      this.cluster_figure.if_load = true;
      // if (typeof this.task_id === 'object') {taskId =undefined}

      return get("get_figure_UMP", {
        params: {
          taskId: this.task_id ,
        }
      })
        .then((res) => {

          // this.cluster_figure.clustered_center = res.clustered_center
          this.cluster_figure.cluster_colors = res.data.cluster_colors;
          this.cluster_figure.clustered_distribution = res.data.clustered_distribution;

          this.cluster_figure.clustered_data = JSON.parse(res.data.clustered_data);
          this.cluster_figure.legend_label = [];
          // data = JSON.parse(data)

          const pointsByLabel = {}; // 用于按标签存储数据点的对象
          const countMap = new Map(); // 用于按标签存储标签出现的次数
          if (this.task_id ){
            this.annotation.currunt_model = "annotation"
            console.log("-----------------------")
            console.log("taskId:",this.task_id )
          }
          else{
            this.annotation.currunt_model = "seurat_clusters"
          }
          console.log("-----------------------")
          console.log("this.annotation.currunt_model:",this.annotation.currunt_model)
            

          // 遍历数据，按标签分组数据点
          this.cluster_figure.clustered_data.forEach((row) => {
            const label = row[this.annotation.currunt_model]; // 根据你的CSV文件实际列名来修改
            const label_go = row["seurat_clusters"];

            if (!pointsByLabel[label]) {
              pointsByLabel[label] = [];
              this.cluster_figure.legend_label.push(label); // 如果标签不存在，则创建一个新的数组并记录顺序
              countMap.set(label, 0); // 统计次数
            }

            // 检查 this.go.clusters 中是否已经存在相同的 label_go
            if (!this.go.clusters.some((cluster) => cluster.value === label_go)) {
              this.go.clusters.push({ label: label_go, value: label_go }); // 添加唯一的 cluster
            }

            countMap.set(label, countMap.get(label) + 1);
            // 将数据点（第一列和第二列的值）添加到对应标签的数组中
            pointsByLabel[label].push([
              parseFloat(row["umap_data.umap_1"]),
              parseFloat(row["umap_data.umap_2"]),
            ]);
          });
          // 对 this.go.clusters 按照 value 进行排序
          this.go.clusters.sort((a, b) => Number(a.value) - Number(b.value)); // 按照数值大小排序

          // 确定比例图数据
          this.scale_diagram.series=res.data.ann_model
          this.fetchScaleDiagramData();
          console.log(
            "this.cluster_figure.legend_label",
            this.cluster_figure.legend_label
          );

          // 计算每个聚类的中心点（简单的平均值）
          this.cluster_figure.clusters = Object.keys(pointsByLabel).map((label) => ({
            label, // 聚类标签
            center: pointsByLabel[label]
              .reduce((acc, [x, y]) => [acc[0] + x, acc[1] + y], [0, 0])
              .map((sum) => sum / pointsByLabel[label].length), // 聚类中心点
            points: pointsByLabel[label], // 聚类中的数据点
          })); // this.cluster_figure.clusters = [label,center,points]

          // 渲染图表
          this.renderChart();
          this.cluster_figure.if_load = false;
        })
        .catch((error) => {
          console.log("error", error);
          if (error.response) {
            ElMessage.error(error.response.data.msg);

            // console.error("qwer", error);
          } else {
            console.error("请求错误:", error.message);
          }
        });
    },

    // 渲染UMAP图
    renderChart() {
      if (!this.$refs.chart) return;
      // 初始化ECharts实例，添加主题
      this.cluster_figure.chart = markRaw(echarts.init(this.$refs.chart));
      console.log("初始化ECharts实例", this.cluster_figure.clusters[0].center[1]);

      // 改进的ECharts配置项
      const option = {
        // backgroundColor: '#f8f9fa',  // 浅灰背景色提升专业感
        title: {
          text: 'UMAP细胞类型分布图',
          left: 'center',
          top: 10,
          textStyle: {
            fontSize: 18,
            fontWeight: 'bold'
          }
        },
        grid: {
          left: '5%',
          right: '5%',
          top: '25%',  // 增加顶部空间，为图例腾出位置
          bottom: '10%',
          containLabel: true
        },
        xAxis: {
          type: "value",
          // 去除坐标轴名称
          name: '',
          // 去除坐标轴线
          axisLine: {
            show: false
          },
          // 去除刻度线
          axisTick: {
            show: false
          },
          // 去除刻度标签
          axisLabel: {
            show: false
          },
          // 保留网格线但使其更淡
          splitLine: {
            show: true,
            lineStyle: {
              type: 'dashed',
              opacity: 0.2
            }
          }
        },
        yAxis: {
          type: "value",
          // 去除坐标轴名称
          name: '',
          // 去除坐标轴线
          axisLine: {
            show: false
          },
          // 去除刻度线
          axisTick: {
            show: false
          },
          // 去除刻度标签
          axisLabel: {
            show: false
          },
          // 保留网格线但使其更淡
          splitLine: {
            show: true,
            lineStyle: {
              type: 'dashed',
              opacity: 0.2
            }
          }
        },
        tooltip: {
          trigger: 'item',
          backgroundColor: 'rgba(255,255,255,0.9)',
          borderColor: '#ccc',
          borderWidth: 1,
          textStyle: {
            color: '#333'
          },
          formatter: function(params) {
            return `${params.seriesName}`;
          }
        },
        dataZoom: [
          {
            type: "inside",
            xAxisIndex: [0],
            start: 0,
            end: 100,
            zoomLock: false,
            filterMode: 'filter'
          },
          {
            type: "inside",
            yAxisIndex: [0],
            start: 0,
            end: 100,
            zoomLock: false,
            filterMode: 'filter'
          },
          {
            type: 'slider',
            xAxisIndex: [0],
            bottom: '1%',
            height: 20,
            start: 0,
            end: 100,
            borderColor: 'transparent',
            fillerColor: 'rgba(80,80,80,0.1)',
            handleStyle: {
              color: '#888'
            },
            // 隐藏滑块上的文字标签
            labelFormatter: function() {
              return '';
            }
          }
        ],
        legend: {
          // 将图例移到顶部
          type: 'scroll',
          orient: 'horizontal',  // 水平布局
          top: 50,  // 位于标题下方
          
          left: 'center',  // 居中显示
          width: '90%',  // 宽度占比
          itemWidth: 15,
          itemHeight: 10,
          textStyle: {
            fontSize: 12
          },
          pageButtonPosition: 'end',  // 翻页按钮位置
          pageButtonGap: 5,  // 翻页按钮与图例的间距
          pageButtonItemGap: 5,  // 翻页按钮之间的间距
          pageIconColor: '#888',  // 翻页按钮颜色
          pageIconInactiveColor: '#ccc',  // 不可用翻页按钮颜色
          pageIconSize: 12,  // 翻页按钮大小
          data: this.cluster_figure.legend_label.sort((a, b) => Number(a) - Number(b)),
          formatter: function(name) {
            // 限制图例文本长度，避免过长
            return name.length > 20 ? name.substring(0, 18) + '...' : name;
          }
        },
        series: this.cluster_figure.clusters.map((cluster, index) => ({
          name: cluster.label,
          type: "scatter",
          data: cluster.points,
          symbolSize: 5,  // 增大点的大小
          symbol: 'circle',
          emphasis: {
            focus: 'series',
            itemStyle: {
              shadowBlur: 10,
              shadowColor: 'rgba(0, 0, 0, 0.1)'
            }
          },
          itemStyle: {
            color: this.cluster_figure.cluster_colors[index],
            opacity: 0.8,  // 添加透明度，减少重叠点的视觉混乱
            borderColor: '#fff',
            borderWidth: 0.3,
            shadowBlur: 2,
            shadowColor: 'rgba(0, 0, 0, 0.2)'
          },
          // markPoint: {
          //   symbol: 'pin',  // 使用大头针样式的标记点
          //   symbolSize: [40, 40],
          //   itemStyle: {
          //     color: this.cluster_figure.cluster_colors[index],
          //     opacity: 0.9,
          //     borderColor: '#fff',
          //     borderWidth: 2,
          //     shadowBlur: 5,
          //     shadowColor: 'rgba(0, 0, 0, 0.3)'
          //   },
          //   data: [
          //     {
          //       name: cluster.label,
          //       coord: cluster.center,
          //       label: {
          //         show: true,
          //         formatter: cluster.label,
          //         position: 'inside',
          //         color: '#fff',  // 白色文字更清晰
          //         fontSize: 13,  // 增大字体大小
          //         fontWeight: 'bold',
          //         textShadowColor: 'rgba(0, 0, 0, 0.6)',  // 增强文字阴影
          //         textShadowBlur: 2,  // 增加文字阴影模糊半径
          //         textShadowOffsetX: 1,  // 添加文字阴影水平偏移
          //         textShadowOffsetY: 1,  // 添加文字阴影垂直偏移
          //       }
          //     }
          //   ]
          // }
        })),


        animation: true,
        animationDuration: 1000,
        animationEasing: 'cubicOut'
      };

      // 渲染图表
      this.cluster_figure.chart.setOption(option, true);
      
      // 添加图表响应式调整
      window.addEventListener('resize', () => {
        this.cluster_figure.chart && this.cluster_figure.chart.resize();
      });
    },
    
    // 获取QC图
    async getQCFigure() {
      try {
        const res = await get("get_figure_QC_final");
        this.imgPaths = res.data.png_paths;
        console.log("this.imgPaths", this.imgPaths);
      } catch (error) {
        if (error.response) {
          console.log("请求失败:", error.response.status);
          console.log("错误响应数据:", error.response.data);
          ElMessage.error(error.response.data.msg || "获取QC终图失败");
        } else {
          console.error("请求错误:", error.message);
        }
        this.imgPaths = [];
      }
    },

    
    // 细胞注释
    // 得到所有的注释方法
    fetchAnnotationModels() {
      this.annotation.value = "";
      get("get_annotation_models")
        .then((response) => {
          // 找到已存在的 celltypist 索引
          const existingIndex = this.annotation.models.findIndex(
            (model) => model.value === "celltypist"
          );

          // 如果存在则替换，不存在则添加
          if (existingIndex !== -1) {
            this.annotation.models[existingIndex] = response.data.celltypist;
          } else {
            this.annotation.models.push(response.data.celltypist);
          }
        })
        .catch((error) => {
          console.error("Error fetching annotation models:", error);
        });
    },
    // 进行注释
    handleSelectChange_annotation(value) {
      this.annotation.value = value;
      this.annotation.currunt_model = "annotation";
      this.annotation.placeholder = this.annotation.value;
      this.cluster_figure.if_load = true;

      if (value === "train") {
        this.annotation.dialogVisible = true; // 如果选择训练，弹出对话框

        this.cluster_figure.if_load = false;
      } else {
        // 处理其他选项
        console.log("选择了:", value);
        const formData = new FormData();
        formData.append("annotation_model", value);
        post("select_annotation_model", formData)
          .then((response) => {
            this.task_id = response.data.task_id_new
            this.getUmapData(); // 更新umap图
          })
          .catch((error) => {
            console.log("error", error);
            if (error.response) {
              console.log("请求失败:", error.response.status);
              console.log("错误响应数据:", error.response.data);
              ElMessage.error(error.response.data.msg);
            } else {
              ElMessage.error("请求错误:", error.message);
              console.error("请求错误:", error.message);
            }
          });
      }
    },

    handleReferenceDataSuccess(response, file) {
      this.annotation.save_path = response.save_path;
      this.$message({
        showClose: true,
        message: "数据上传成功",
        type: "success",
      });
    },
    async submitReferenceData() {
      try {
        this.annotation.dialogVisible = false;
        ElNotification({
          title: "Info",
          message: "数据已提交，等训练完成可以在下拉框中找到该模型",
          type: "info",
        });
        const formData = new FormData();
        formData.append("save_path", this.annotation.save_path);
        formData.append("model_name", this.annotation.modelName);
        const response = await post("submit_reference_data", formData);
      } catch (error) {
        // 优化错误提示
        const errorMsg = error.response?.data?.message || error.message || "未知错误";
        this.$message.error("提交失败：" + errorMsg);
      } finally {
        // 确保在任何情况下都关闭loading
        console.log("");
      }
    },

    async LargeModelAnnotation() {
      if (!this.task_id) {
        this.$message({
          showClose: true,
          message: "还没有选择任务，请于项目管理页面选择任务",
          type: "error",
        });
        this.$router.push("/upload");
        return;
      }
      try {
        this.annotation.LargeModelDialogVisible = true

        const response = await post("large_model_annotation", {
          taskId: this.task_id ,
        });
        this.annotation.tableData = response.data.anno;
      } catch (error) {
        // 检查后端返回的错误信息
        const msg = error.response?.data?.message || error.message || "未知错误";
        this.$message({
          showClose: true,
          message: msg,
          type: "error",
        });

      }
    },

    // 计算marker基因
    async CalculateMarkerGene() {
      this.marker.if_load = true;
      try {
        await post("calculate_marker_gene", {});
      } catch (error) {
        console.error("计算marker基因失败:", error);
        ElMessage.error("计算失败");
      } finally {
        this.marker.if_load = false;
      }
    },

    // 高可变基因
    // 高可变基因数据获取
    handleSelectChange_marker(value) {
      const formData = new FormData();
      formData.append("num", value);
      formData.append("taskId", this.task_id);

      post("get_marker_genes", formData)
        .then((response) => {
          this.marker.tableData = response.data.gene_list;
        })
        .catch((error) => {
          console.error("获取marker基因失败:", error);
          ElMessage.error("获取失败");
        });
    },

    // 下载表格数据为CSV文件
    exportToCSV() {
      // 检查 tableData 是否存在且为数组
      if (!Array.isArray(this.marker.tableData)) {
        console.error("tableData 不是一个数组");
        return;
      }

      let csvContent = "gene,p_val,avg_log2FC,pct.1,pct.2,p_val_adj,cluster\n";
      this.marker.tableData.forEach((row) => {
        // 确保每一行都有必要的属性
        csvContent += `${row.gene || "NA"},${row.p_val || "NA"},${
          row.avg_log2FC || "NA"
        },${row["pct.1"] || "NA"},${row["pct.2"] || "NA"},${row.p_val_adj || "NA"},${
          row.cluster || "NA"
        }\n`;
      });

      const blob = new Blob([csvContent], { type: "text/csv;charset=utf-8;" });
      const link = document.createElement("a");
      link.href = URL.createObjectURL(blob);
      link.download = "table-data.csv";
      // 触发下载
      document.body.appendChild(link); // 有时需要先将元素添加到 DOM 中才能触发点击事件
      link.click();
      // 清理：移除链接并从 DOM 中删除
      document.body.removeChild(link);
      URL.revokeObjectURL(link.href);
    },

    // 绘制DotPlot图
    async getDotPlotFigure(value) {
      const formData = new FormData();
      formData.append("num", value);
      formData.append("taskId", this.task_id);

      try {
        const response = await post("get_figure_dotplot", formData);
        this.marker.dotplot_data = JSON.parse(response.data.dotplot_data);
        this.initDotPlotChart();
      } catch (error) {
        console.error("获取气泡图失败:", error);
        ElMessage.error("获取失败");
      }
    },
    initDotPlotChart() {
      // 初始化两个图表实例
      const xAxisChart = echarts.init(this.$refs.xAxisChart);
      const chart_dotplot = echarts.init(this.$refs.chart_dotplot);

      const annotations = [
        ...new Set(this.marker.dotplot_data.map((item) => item.Annotation)),
      ];
      const genes = [...new Set(this.marker.dotplot_data.map((item) => item.Gene))];

      const seriesData = this.marker.dotplot_data.map((item) => ({
        name: item.Gene,
        value: [
          annotations.indexOf(item.Annotation),
          genes.indexOf(item.Gene),
          item.Mean_Expression,
          item.Expression_Percentage,
        ],
      }));

      // X轴配置
      const xAxisOption = {
        grid: {
          // top: 60,
          left: 0,
          bottom: 1,
          // left: 0
        },
        xAxis: {
          type: "category",
          data: annotations,
          position: "top",
          axisLabel: {
            rotate: 45,
            interval: 0,
            show: true,
            margin: -10,
          },
          axisTick: {
            alignWithLabel: true,
            show: true,
          },
          axisLine: {
            show: true,
            lineStyle: {
              color: "#333",
            },
          },
          z: 3,
        },
        yAxis: {
          show: false,
        },
        series: [],
      };

      // 合并后的图表配置
      const combinedOption = {
        grid: {
          top: 0,
          right: 100,
          bottom: 0,
          left: 150,
        },
        xAxis: {
          type: "category",
          data: annotations,
          show: false,
        },
        yAxis: {
          type: "category",
          data: genes,
          inverse: true,
          position: "left", // Y 轴位于左侧
          axisLabel: {
            interval: 0,
            show: true,
            align: "right",
          },
          axisTick: {
            alignWithLabel: true,
            show: true,
          },
          axisLine: {
            show: true,
          },
        },
        tooltip: {
          trigger: "item",
          backgroundColor: "rgba(255,255,255,0.95)",
          borderColor: "#ccc",
          borderWidth: 1,
          padding: 10,
          textStyle: {
            color: "#333",
            fontSize: 14,
          },
          formatter: function (params) {
            const annotation = annotations[params.value[0]];
            const gene = genes[params.value[1]];
            const meanExpression = params.value[2].toFixed(3);
            const expressionPercentage = params.value[3].toFixed(3);

            return `
              <div style="text-align: left;">
                <b>Annotation:</b> ${annotation}<br>
                <b>Gene:</b> ${gene}<br>
                <b>Avg Exp:</b> ${meanExpression}<br>
                <b>Pct1:</b> ${expressionPercentage}
              </div>
            `;
          },
        },
        visualMap: [
          {
            // 平均表达值映射到颜色
            type: "continuous",
            dimension: 2, // value[2] 是 Mean_Expression
            min: 0,
            max: 1,
            calculable: true,
            orient: "vertical",
            left: "right",
            top: "20%",
            inRange: {
              color: ["#313695", "#ffffff", "#d73027"], // 蓝到红渐变
            },
            text: ["高表达", "低表达"],
            textStyle: {
              color: "#333",
            },
          },
          {
            // 表达百分比映射到大小
            type: "continuous",
            dimension: 3, // value[3] 是 Expression_Percentage
            min: 0,
            max: 0.5,
            calculable: true,
            orient: "vertical",
            left: "right",
            top: "top",
            inRange: {
              symbolSize: [2, 20], // 圆点大小范围
            },
            text: ["高比例", "低比例"],
            textStyle: {
              color: "#333",
            },
          },
        ],
        series: [
          {
            type: "scatter",
            data: seriesData,
            itemStyle: {
              opacity: 0.8,
            },
          },
        ],
      };

      // 设置图表的配置
      xAxisChart.setOption(xAxisOption);
      chart_dotplot.setOption(combinedOption);
    },

    // 特征图
    // 得到选择列表
    fetchMarkerGenes() {
      get("get_gene_list")
        .then((response) => {
          this.marker.gene_list = response.data.gene_list;
          this.tcga.gene_list = response.data.gene_list_all;
          // this.pseudotime.gene_list = response.data.gene_list;
        })
        .catch((error) => {
          console.error("Error fetching annotation models:", error);
        });
      console.log("this.marker.gene_list", this.marker.gene_list);
    },

    async submitSelection_Significance() {
      if (this.SignificanceQC.checkedCells.length !== 2) {
        this.$message.warning("请选择两个细胞类型进行比较");
        return;
      }

      this.SignificanceQC.loading = true;
      try {
        const response = await axios.post("/api/get_SignificanceQC", {
          selectedCells: this.SignificanceQC.checkedCells,
          selectedGene: this.marker.gene,
        });

        // 使用Vue的响应式更新方式
        this.SignificanceQC.imagePath = response.data.imagePath;
        console.log(this.SignificanceQC.imagePath);
        // this.$set(this.SignificanceQC, 'imagePath', response.data.imagePath);
        this.$message.success("分析完成");
      } catch (error) {
        this.$message.error("提交失败：" + error.message);
      } finally {
        this.SignificanceQC.loading = false;
      }
    },


    // 获取分类模式
    fetchCategoryMode() {
      get("get_category_model")
        .then((response) => {
          this.scale_diagram.models = response.data.model_list;
        })
        .catch((error) => {
          console.error("Error fetching scale_diagram models:", error);
        });
    },
    // 应用分类模式
    applyCategoryMode() {
      if (!this.scale_diagram.categoryMode) {
        this.$message.warning("请输入分类模式");
        return;
      }

        this.fetchScaleDiagramData();
        if (this.scale_diagram.rawData.length > 0) {
        } else {
          this.scale_diagram.myChart.showLoading();
        }
        this.$message.success(`已应用"${this.scale_diagram.categoryMode}"分类模式`);
    },
    // 渲染细胞比例scale_diagram图
    async fetchScaleDiagramData() {
      this.scale_diagram.myChart = echarts.init(this.$refs.chart_scale_diagram);
      this.scale_diagram.myChart.showLoading(); // 开始加载动画
      try {
        
        const response = await get("get_scale_diagram_data", {
          params: {
            categoryMode: this.scale_diagram.categoryMode,
            seriesNames: this.scale_diagram.series,
          },
        });

        // 检查响应数据
        if (response.data.status === "success") {
          this.scale_diagram.rawData = response.data.rawData;
          this.scale_diagram.categories = response.data.categories;
          this.scale_diagram.seriesNames = response.data.seriesNames;
          this.calculateTotalData();
          this.updateChartOption();
        } else {
          this.$message.error(response.data.message || "获取数据失败");
        }
      } catch (error) {
        console.error("获取数据失败:", error);
        this.$message.error("获取比例图数据失败，请稍后重试");
      } finally {
        this.scale_diagram.myChart.hideLoading();
      }
    },
    // 渲染细胞比例scale_diagram图

    // 计算每个类别的总数据
    calculateTotalData() {
      this.scale_diagram.totalData = [];
      if (
        this.scale_diagram.rawData.length > 0 &&
        this.scale_diagram.rawData[0].length > 0
      ) {
        for (let i = 0; i < this.scale_diagram.rawData[0].length; ++i) {
          let sum = 0;
          for (let j = 0; j < this.scale_diagram.rawData.length; ++j) {
            sum += this.scale_diagram.rawData[j][i];
          }
          this.scale_diagram.totalData.push(sum);
        }
      }
    },

    // 更新图表选项
    updateChartOption() {
      if (this.scale_diagram.rawData.length === 0) return;

      this.scale_diagram.myChart.hideLoading();

      const maxLen = Math.max(...this.scale_diagram.seriesNames.map(n => n.length));
      const legendWidth = maxLen * 7; // 7px/字符，40px 图标+留白
      const grid = {
        left: 100,
        right: 120, // 再加 20px 安全间隙
        top: 50,
        bottom: 50,
      };


      let series;
      if (this.scale_diagram.selectedView === "percentage") {
        series = this.scale_diagram.seriesNames.map((name, sid) => ({
          name,
          type: "bar",
          stack: "total",
          barWidth: "80%",
          emphasis: {
            focus: "series",
          },
          label: {
            show: true,
            formatter: (params) => Math.round(params.value * 1000) / 10 + "%",
          },
          data: this.scale_diagram.rawData[sid].map((d, did) =>
            this.scale_diagram.totalData[did] <= 0
              ? 0
              : (d / this.scale_diagram.totalData[did]).toFixed(3)
          ),
          itemStyle: {
            color: this.cluster_figure.cluster_colors[sid],
          },
        }));
      } else {
        series = this.scale_diagram.seriesNames.map((name, sid) => ({
          name,
          type: "bar",
          stack: "total",
          barWidth: "80%",
          emphasis: {
            focus: "series",
          },
          label: {
            show: true,
            formatter: (params) => params.value,
          },
          data: this.scale_diagram.rawData[sid],
        }));
      }

      const option = {
        title: {
          text:
            this.scale_diagram.selectedView === "percentage"
              ? "各类别占比分布"
              : "各类别绝对值分布",
          left: "center",
        },
        tooltip: {
          trigger: "axis",
          axisPointer: {
            type: "shadow",
          },
        },
        legend: {
          type: "scroll",
          orient: "vertical",
          right: 10,
          top: "middle", // 垂直居中
          itemGap: 5, // 图例项间距
          itemWidth: 14, // 图例标记宽度
          itemHeight: 14, // 图例标记高度
          selectedMode: false, // 添加这行来禁用图例的点击交互
          textStyle: {
            fontSize: 12,
            rich: {
              // 设置文本超出省略
              ellipsis: {
                width: 100, // 最大宽度
                overflow: "truncate", // 超出部分省略
                ellipsis: "...",
              },
            },
          },
          formatter: function (name) {
            // 限制图例名称长度，超过15个字符显示省略号
            return name.length > 15 ? name.substring(0, 12) + "..." : name;
          },
        },
        grid,
        yAxis: {
          type: "value",
          max: this.scale_diagram.selectedView === "percentage" ? 1 : undefined, // 添加这行，百分比视图时设置最大值为1
          axisLabel: {
            formatter:
              this.scale_diagram.selectedView === "percentage"
                ? (value) => (value * 100).toFixed(0) + "%"
                : (value) => value,
          },
        },
        xAxis: {
          type: "category",
          data: this.scale_diagram.categories,
        },
        series,
      };

      this.scale_diagram.myChart.setOption(option, true); // 添加 true 参数以完全替换配置
    },
    // 更改类别视图
    handleViewChange() {
      this.updateChartOption();
    },

    // 渲染Go富集条形图
    async handleSelectChange_go() {
      this.$nextTick(() => {
        if (!this.$refs.chart_bar) {
          console.error('图表容器不存在');
          return;
        }

        // 如果已存在实例，先销毁
        if (this.go.chart) {
          this.go.chart.dispose();
        }

        // 初始化图表实例
        this.go.chart = echarts.init(this.$refs.chart_bar);
        this.go.chart.showLoading();

        const formData = new FormData();
        formData.append("choose_cluster", this.go.choose_cluster);

        post("get_figure_barplot", formData)
          .then((response) => {
            this.go.bar_data = response.data.bar_data;
            this.initBarChart();
          })
          .catch((error) => {
            console.error("获取GO图失败:", error);
            ElMessage.error("获取失败");
          })
          .finally(() => {
            if (this.go.chart) {
              this.go.chart.hideLoading();
            }
          });
      });
    },
    async fetchGoFigure() {
      this.$nextTick(() => {
        if (!this.$refs.chart_bar) {
          console.error('图表容器不存在');
          return;
        }

        // 如果已存在实例，先销毁
        if (this.go.chart) {
          this.go.chart.dispose();
        }

        // 初始化图表实例
        this.go.chart = echarts.init(this.$refs.chart_bar);
        this.go.chart.showLoading();

        const formData = new FormData();
        formData.append("taskId", this.task_id);

        post("fetch_go_figure", formData)
          .then((response) => {
            this.go.bar_data = response.data.bar_data;
            this.initBarChart();
          })
          .catch((error) => {
            console.error("获取GO图失败:", error);
            ElMessage.error("获取失败");
          })
          .finally(() => {
            if (this.go.chart) {
              this.go.chart.hideLoading();
            }
          });
      });
    },
    // 初始化两个图表实例
    initBarChart() {
      var option;
      // 处理 bar_data 以适应 ECharts 的数据格式
      const data = this.go.bar_data.map((item) => {
        return {
          value: item[0], // Count
          name: item[1], // Description
          pAdjust: item[2], // p.adjust
        };
      });

      // 获取 p.adjust 的最大值和最小值
      const pAdjustValues = data.map((item) => item.pAdjust);
      const minPAdjust = Math.min(...pAdjustValues);
      const maxPAdjust = Math.max(...pAdjustValues);

      option = {
        tooltip: {
          trigger: "item", // 触发类型为 item
          formatter: (params) => {
            // 显示 Description 和 p.adjust，以科学计数法格式化
            return `${params.data[1]}: ${params.data[2].toExponential(2)}`; // params.data[1] 是 Description，params.data[2] 是 p.adjust
          },
        },
        dataset: {
          source: data.map((item) => [item.value, item.name, item.pAdjust]), // 也保留 p.adjust
        },
        grid: {
          containLabel: true,
          top: "5%", // 增大上边距
          bottom: "5%", // 增大下边距
          right: 90, // 增大下边距
          left: "5%", // 增大下边距
        },
        
        xAxis: {
          name: "Count",
          nameLocation: "center",
          nameGap: 20,
        },
        yAxis: {
          type: "category",
          inverse: true, // 让数据从上往下排列
          // axisLabel: {
          //   formatter: (value) => {
          //     // 设置最大长度
          //     const maxLength = 10;
          //     // 如果长度超过最大长度，则截断并添加省略号
          //     if (value.length > maxLength) {
          //       return value.substring(0, maxLength) + "..."; // 截断并添加省略号
          //     }
          //     return value; // 否则返回原值
          //   },
          // },
          axisLabel: {
            /*  ***删除原来的 formatter***  */
            width: 150,          // 给标签预留总宽度（px）
            overflow: 'break',   // echarts 5.4+ 自动换行
            interval: 0,         // 强制显示全部标签
            fontSize: 12,
            lineHeight: 10,      // 行高，保证两行美观
            rich: {              // 统一两行样式
              desc: {
                fontSize: 12,
                lineHeight: 10
              }
            }
          }
        },
        visualMap: {
          show: true,
          dimension: 2, // 使用 p.adjust 作为映射维度
          min: minPAdjust, // 设置最小值
          max: maxPAdjust, // 设置最大值
          inRange: {
            color: ["#df6664", "#a2758f", "#367eb9"],
          },
          outOfRange: {
            color: "#ccc", // 超出范围的颜色
          },
          text: [
            `Max: ${maxPAdjust.toExponential(2)}`,
            `Min: ${minPAdjust.toExponential(2)}`,
          ], // 显示具体的最小值和最大值
          // calculable: true, // 允许用户调整范围
          right: 0,
          top: "center",
          itemHeight: 100, // 设置图例的高度
          itemWidth: 20, // 设置图例的宽度
          inverse: true, // 反转颜色映射
        },
        series: [
          {
            type: "bar",
            encode: {
              x: "Count",
              y: "Description",
            },
            barWidth: 28, // 设置条的宽度，可以根据需要调整
          },
        ],
      };

      option && this.go.chart.setOption(option);
    },

    
    // 工具函数
    // downloadImage() {
    //   if (!this.infercnv.imgPath) return;

    //   const imgUrl = require("../assets/result/" + this.infercnv.imgPath);
    //   const suffix = this.infercnv.imgPath.slice(this.infercnv.imgPath.lastIndexOf("."));
    //   const filename = "infercnv_" + Date.now() + suffix;

    //   fetch(imgUrl)
    //     .then((response) => response.blob())
    //     .then((blob) => {
    //       const blobUrl = URL.createObjectURL(new Blob([blob]));
    //       const link = document.createElement("a");
    //       link.href = blobUrl;
    //       link.download = filename;
    //       document.body.appendChild(link);
    //       link.click();
    //       URL.revokeObjectURL(blobUrl);
    //       link.remove();
    //     })
    //     .catch((error) => {
    //       this.$message.error("下载失败：" + error.message);
    //     });
    // },


    submitTCGAForm(formName) {
      if (!this.task_id) {
        this.$message({
          showClose: true,
          message: "还没有选择需要重绘的任务，请于项目管理页面选择任务",
          type: "error",
        });
        this.tcga.dialogFormVisible = false;
        this.$router.push("/upload");
        return;
      }
      //创建formdata对象
      const formData = new FormData();
      // 添加用户ID到表单数据
      formData.append("task_id", this.task_id);
      formData.append("outlier_size", this.tcga.form.outlier_size);
      formData.append("frame_size", this.tcga.form.frame_size);
      formData.append("box_width", this.tcga.form.box_width);
      formData.append("axis_textsize", this.tcga.form.axis_textsize);
      formData.append("base_size", this.tcga.form.base_size);
      formData.append("ylimits", this.tcga.form.ylimits);

      this.$refs[formName].validate((valid) => {
        if (valid) {
          this.$message("正在重绘中，请稍等……");
          this.tcga.dialogFormVisible = false;

          // 使用 api.js 中的 post 方法
          post("tcga_resize", formData)
            .then((response) => {
              console.log(response.data);
              this.QueryLatestResult_tcga()
            })
            .catch((error) => {
              this.$message({
                showClose: true,
                message: "重绘失败：" + error.message,
                type: "error",
              });
            });
        } else {
          console.log("重绘失败");
          return false;
        }
      });
    },

    downloadImage(url) {
      if (!url) return;
      console.log(url);

      // 从URL中提取文件名
      const filename = url.split("/").pop() || "image_" + Date.now() + ".jpeg";

      // 使用原始URL直接下载
      fetch(url)
        .then((response) => {
          if (!response.ok) {
            throw new Error("Network response was not ok");
          }
          return response.blob();
        })
        .then((blob) => {
          const blobUrl = window.URL.createObjectURL(blob);
          const link = document.createElement("a");
          link.href = blobUrl;
          link.download = filename;
          document.body.appendChild(link);
          link.click();
          document.body.removeChild(link);
          window.URL.revokeObjectURL(blobUrl);
        })
        .catch((error) => {
          console.error("Download failed:", error);
          ElMessage.error("下载失败：" + error.message);
        });
    },
    // 添加跳转方法
    jumpToTutorial(section) {
      // 修改路径为根路径
      window.open(`/tutorial.html#${section}`, "_blank");
    },
  },
  mounted() {
    this.getQCFigure();
    this.fetchCategoryMode();
    this.fetchMarkerGenes();
    this.getCellList();
  },

  beforeUnmount() {
    if (this.cellchat.myChart) {
      this.cellchat.myChart.dispose();
      this.cellchat.myChart = null;
    }
    if (this.scale_diagram.myChart) {
      this.scale_diagram.myChart.dispose();
      this.scale_diagram.myChart = null;
    }
    if (this.go.chart) {
        this.go.chart.dispose();
        this.go.chart = null;
      }
  },
};
</script>

<style scoped>
.insert {
  width: 30%;
}

.demo-image__lazy {
  display: block;
  min-height: 200px;
  height: 100%;
  justify-content: space-between;
  margin-top: 0px;
  display: flex;
}

.image-slot {
  display: flex;
  justify-content: center;
  align-items: center;
  width: 100%;
  height: 100%;
  background: var(--el-fill-color-light);
  color: var(--el-text-color-secondary);
  font-size: 30px;
}

.znjy-container {
  padding: 0px;
  display: flex;
  flex-direction: column;
  align-items: center;
}

.chart-container {
  margin: 5px 0;
  box-shadow: 0 2px 12px 0 rgba(0, 0, 0, 0.1);
  border-radius: 4px;
  padding: 5px;
}

.category-input {
  display: flex;
  gap: 10px;
}

.category-input-field {
  width: 100px;
}

.controls {
  display: flex;
  justify-content: flex-start;
  /* 从左到右排列 */
  align-items: center;
  /* 垂直居中 */
  gap: 5px;
  /* 控制元素之间的间距 */
  margin-top: 10px;
  flex-wrap: nowrap;
  /* 不换行 */
}

.transfer-container {
  display: flex;
  flex-direction: column;
  gap: 20px;
}

.secondary-selection {
  margin-top: 20px;
  padding: 20px;
  border: 1px solid #dcdfe6;
  border-radius: 4px;
}

.el-checkbox-group {
  display: flex;
  flex-wrap: wrap;
  gap: 10px;
}

.result-display {
  margin-top: 20px;
  padding-top: 20px;
  border-top: 1px solid #ebeef5;
}

h3,
h4 {
  margin: 0 0 15px 0;
  color: #303133;
}

/* 添加以下样式 */
:deep(.n-split) {
  width: 100%;
}

.demo-image__placeholder .block {
  padding: 30px 0;
  text-align: center;
  border-right: solid 1px var(--el-border-color);
  display: inline-block;
  width: 49%;
  box-sizing: border-box;
  vertical-align: top;
}

.demo-image__placeholder .demonstration {
  display: block;
  color: var(--el-text-color-secondary);
  font-size: 14px;
  margin-bottom: 20px;
}

.demo-image__placeholder .el-image {
  padding: 0 5px;
  max-width: 300px;
  max-height: 200px;
}

.demo-image__placeholder.image-slot {
  display: flex;
  justify-content: center;
  align-items: center;
  width: 100%;
  height: 100%;
  background: var(--el-fill-color-light);
  color: var(--el-text-color-secondary);
  font-size: 14px;
}

.demo-image__placeholder .dot {
  animation: dot 2s infinite steps(3, start);
  overflow: hidden;
}

.demo-image__lazy {
  height: 400px;
  overflow-y: auto;
}

.demo-image__lazy .el-image {
  display: block;
  min-height: 200px;
  margin-bottom: 10px;
}

.demo-image__lazy .el-image:last-child {
  margin-bottom: 0;
}

.photo-gallery {
  padding: 20px;
  width: 100%;
  max-width: 1200px;
  margin: 0 auto;
}

.gallery-row {
  display: flex;
  justify-content: space-between;
  margin-bottom: 30px;
  gap: 20px;
}

.photo-item {
  flex: 1;
  text-align: center;
}

.photo-item h3 {
  margin-bottom: 10px;
  color: #333;
}

.el-image {
  width: 100%;
  height: 300px;
  border-radius: 8px;
  box-shadow: 0 2px 12px 0 rgba(0, 0, 0, 0.1);
}

.image-slot {
  display: flex;
  justify-content: center;
  align-items: center;
  width: 100%;
  height: 100%;
  background: var(--el-fill-color-light);
  color: var(--el-text-color-secondary);
  font-size: 30px;
}

/* ... 其他样式 ... */
.empty-placeholder {
  display: flex;
  flex-direction: column;
  justify-content: center;
  align-items: center;
  height: 200px;
  color: var(--el-text-color-secondary);
}

.empty-placeholder .el-icon {
  margin-bottom: 10px;
}

.empty-placeholder p {
  margin: 0;
  font-size: 14px;
}

.dialog-footer {
  display: flex;
  justify-content: flex-end;
}
</style>
