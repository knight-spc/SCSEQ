import axios from 'axios'
import { ElMessage } from 'element-plus'

const instance = axios.create({
  baseURL: '/api',
  timeout: 7200000  // 设置为2小时 (2 * 60 * 60 * 1000 毫秒)
  })

// 添加请求拦截器
instance.interceptors.request.use(
  config => {
    const userId = localStorage.getItem('user_id')
    const itemId = localStorage.getItem('item_id')
    if (!config.url.includes('login')) {
      if (!userId) {
        ElMessage.warning('请先登录')
        return Promise.reject('未登录')
      }
      config.headers['User-Id'] = userId

      if (itemId) {
        config.headers['Item-Id'] = itemId
      }
    }
    return config
  },
  error => {
    return Promise.reject(error)
  }
)

export function  get(url, params = {}) {
    return instance.get(url, {
      params
    })
  }

export function  post(url, data = {}) {
    return instance.post(url, data)
  }

// 登陆
export function reqLogin(username, password){
  return post('login',{username,password})
}

export function get_datastream(data_name, start_time, end_time, limit){
  ElMessage({
    message: '获取数据成功',
    type: 'success',
  })
  return get('get_datastream',{data_name, start_time, end_time, limit})
}

// export function get_file(filename){
//   return axios.get('api/get_file', {
//     params: {
//       filename:filename
//     },
//     responseType: 'arraybuffer' // 将响应类型设置为arraybuffer
//   })
// }

// export function generate_excel(data_name, start_time, end_time, limit){
//   return axios.get('api/generate_excel', {
//     params: {
//       data_name:data_name,
//       start_time:start_time,
//       end_time:end_time,
//       limit:limit
//     },
//     responseType: 'arraybuffer' // 将响应类型设置为arraybuffer
//   })
// }

export function send_move(message){
  return get(`send_move?message=${message}`)
}