<template>
  <div class="hero">
  <h1 class="title">

       <span>单细胞转录组</span>
       <span>测序分析云平台</span>
     </h1>
   </div>
     <el-carousel :interval="5000" arrow="always" class="photo">
       <el-carousel-item v-for="item in 2" :key="item">
         <h3 text="2xl" justify="center"></h3>
       </el-carousel-item>
     </el-carousel>
 
</template>

<style scoped>

@import '../assets/font/font.css';
.hero {
/* height: 100px; */
/* background-image: url('../assets/images/background.jpg');   */
background-size: cover;
display: flex;
align-items: center;
justify-content: center;
text-align: center;
}

.title {
color: black;
font-size: 48px;
line-height: 56px;
font-weight: bold;
margin-top:10px;
}

.el-carousel__item h3 {
 /* height: 500px; */

 color: #475669;
 opacity: 0.75;
 line-height: 500px;
 margin: 0;
 text-align: center;
}

.el-carousel__item:nth-child(2n) {

 background-image: url('../assets/images/scseq1.png');
 height: 400px;
 background-repeat: no-repeat;
 background-position: center center;
 background-size: cover;
/* filter: blur(5px); */
}

.el-carousel__item:nth-child(2n + 1) {
 height: 400px;

 background-image: url('../assets/images/scseq2.png');
 background-repeat: no-repeat;
 background-position: center center;
 background-size: cover;
/* filter: blur(5px); */
}

.photo{
 height: 500px;
}
</style>
