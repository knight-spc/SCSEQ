const { defineConfig } = require('@vue/cli-service')
const path = require('path')

module.exports = defineConfig({
  transpileDependencies: true,
  lintOnSave: false,
  productionSourceMap: false, 
  publicPath: './',
  
  devServer: {

    port: 8080, // 端口
    proxy: {
        '/api': { 
            target: '', 
            changeOrigin: true,
            ws: true,
        }
    },
    

  },

  configureWebpack: {
    performance: {
      hints: false
    },
    optimization: {
      splitChunks: {
        chunks: 'all' 
      }
    }
  }
})



  