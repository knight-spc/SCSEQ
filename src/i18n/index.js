import { createI18n } from 'vue-i18n'
import en from './en.json'
import zh from './zh.json'

const i18n = createI18n({
  legacy: false,
  locale: 'en',        // 👈 默认英文（对审稿人非常重要）
  fallbackLocale: 'en',
  messages: { en, zh }
})

export default i18n
