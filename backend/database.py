import pymysql
import json  # 添加 json 模块导入



host = ''


def create_database():
    conn = pymysql.connect(host=host, port=3306, user='sunpc', password='spc20020629', database='scseq')
    c = conn.cursor()

    c.execute("""
        CREATE TABLE IF NOT EXISTS users (
            id INTEGER PRIMARY KEY AUTO_INCREMENT,
            name VARCHAR(50) UNIQUE NOT NULL,
            password VARCHAR(100) NOT NULL,
            phone VARCHAR(20) UNIQUE
        ) default charset=utf8
    """)
    print("Table 'users' checked/created.")

    c.execute("""
        CREATE TABLE IF NOT EXISTS items (
            item_id INTEGER PRIMARY KEY AUTO_INCREMENT,
            user_id INTEGER NOT NULL,
            item_name TEXT,
            save_path TEXT,
            work_path TEXT,
            species TEXT,
            notes TEXT,
            created_time TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
            FOREIGN KEY (user_id) REFERENCES users(id)
        ) default charset=utf8
    """)
    print("Table 'items' checked/created.")

    c.execute("""
        CREATE TABLE IF NOT EXISTS tasks (
            task_id INTEGER PRIMARY KEY AUTO_INCREMENT,
            item_id INTEGER NOT NULL,
            type VARCHAR(50),
            parameter TEXT,
            result_path TEXT,
            time TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
            sh_id INTEGER,  
            FOREIGN KEY (item_id) REFERENCES items(item_id)
        ) default charset=utf8
    """)
    print("Table 'tasks' checked/created.")

    conn.commit()
    conn.close()


def connect_db():
    try:
        conn = pymysql.connect(
            host=host, 
            port=3306, 
            user='sunpc', 
            password='spc20020629', 
            database='scseq'
        )
        return conn
    except Exception as e:
        print(f"数据库连接失败: {str(e)}")
        raise e

def add_item(user_id, item_name, save_path, work_path, species, notes):
    conn = connect_db()
    try:
        with conn.cursor() as cursor:
            sql = """INSERT INTO items 
                    (user_id, item_name, save_path, work_path, species, notes) 
                    VALUES (%s, %s, %s, %s, %s, %s)"""
            cursor.execute(sql, (user_id, item_name, save_path, work_path, species, notes))
        conn.commit()
        return cursor.lastrowid
    finally:
        conn.close()

def get_item(item_id):
    conn = connect_db()
    try:
        with conn.cursor() as cursor:
            sql = "SELECT * FROM items WHERE item_id = %s"
            cursor.execute(sql, (item_id,))
            return cursor.fetchone()
    finally:
        conn.close()
def get_item_workpath(item_id):
    conn = connect_db()
    try:
        with conn.cursor() as cursor:
            sql = """SELECT work_path  
                        FROM items WHERE item_id = %s 
                        ORDER BY item_id DESC"""
            cursor.execute(sql, (item_id,))
            return cursor.fetchone()
    finally:
        conn.close()
def get_item_notes(item_id):
    conn = connect_db()
    try:
        with conn.cursor() as cursor:
            sql = """SELECT notes  
                        FROM items WHERE item_id = %s 
                        ORDER BY item_id DESC"""
            cursor.execute(sql, (item_id,))
            return cursor.fetchone()
    finally:
        conn.close()

def get_user_items(user_id):
    conn = connect_db()
    try:
        with conn.cursor() as cursor:
            sql = "SELECT * FROM items WHERE user_id = %s"
            cursor.execute(sql, (user_id,))
            return cursor.fetchall()
    finally:
        conn.close()

def update_item(item_id, save_path=None, parameter=None, result_img=None, calculation_results=None):
    conn = connect_db()
    try:
        with conn.cursor() as cursor:
            updates = []
            params = []
            if save_path:
                updates.append("save_path = %s")
                params.append(save_path)
            if parameter:
                updates.append("parameter = %s")
                params.append(parameter)
            if result_img:
                updates.append("result_img = %s")
                params.append(result_img)
            if calculation_results:
                updates.append("calculation_results = %s")
                params.append(calculation_results)
            if updates:
                sql = f"UPDATE items SET {', '.join(updates)} WHERE item_id = %s"
                params.append(item_id)
                cursor.execute(sql, tuple(params))
                conn.commit()
                return True
        return False
    finally:
        conn.close()

def delete_item(item_id):
    conn = connect_db()
    try:
        with conn.cursor() as cursor:
            sql = "DELETE FROM items WHERE item_id = %s"
            cursor.execute(sql, (item_id,))
        conn.commit()
        return True
    except Exception as e:
        print(f"删除任务失败: {e}")
        return False
    finally:
        conn.close()

def get_user_tasks(user_id):
    conn = connect_db()
    try:
        with conn.cursor() as cursor:
            sql = "SELECT * FROM tasks WHERE user_id = %s"
            cursor.execute(sql, (user_id,))
            return cursor.fetchall()
    finally:
        conn.close()

def update_item_sh_id(job_id, job_id_new):
    conn = connect_db()
    try:
        with conn.cursor() as cursor:
            sql = "UPDATE tasks SET sh_id = %s WHERE sh_id = %s"
            cursor.execute(sql, (job_id_new, job_id))
        conn.commit()
        return True
    except Exception as e:
        print(f"更新sh_id失败: {e}")
        return False
    finally:
        conn.close()
def delete_task(task_id):
    conn = connect_db()
    try:
        with conn.cursor() as cursor:
            sql = "DELETE FROM tasks WHERE task_id = %s"
            cursor.execute(sql, (task_id,))
        conn.commit()
        return True
    except Exception as e:
        print(f"删除任务失败: {e}")
        return False
    finally:
        conn.close()

def add_new_task(item_id, task_type, parameter=None, result_path=None, sh_id=None):
    conn = connect_db()
    try:
        if isinstance(parameter, dict):
            parameter = json.dumps(parameter)
        if isinstance(result_path, dict):
            result_path = json.dumps(result_path)
        with conn.cursor() as cursor:
            sql = """INSERT INTO tasks 
                    (item_id, type, parameter, result_path, sh_id) 
                    VALUES (%s, %s, %s, %s, %s)"""
            cursor.execute(sql, (item_id, task_type, parameter, result_path, sh_id))
        conn.commit()
        return cursor.lastrowid  
    except Exception as e:
        print(f"添加任务失败: {e}")
        return None
    finally:
        conn.close()

def get_latest_task_by_type(item_id, task_type):

    conn = connect_db()
    try:
        with conn.cursor() as cursor:
            sql = """SELECT * FROM tasks 
                    WHERE item_id = %s AND type = %s 
                    ORDER BY time DESC 
                    LIMIT 1"""
            cursor.execute(sql, (item_id, task_type))
            result = cursor.fetchone()
            if result:
                result = list(result)
                try:
                    if isinstance(result[3], str):
                        result[3] = json.loads(result[3])
                    if isinstance(result[4], str):
                        result[4] = json.loads(result[4])
                except json.JSONDecodeError:
                    pass  
                return tuple(result)
            return None
    finally:
        conn.close()

def get_task_by_id(task_id):
    conn = connect_db()
    try:
        with conn.cursor() as cursor:
            sql = "SELECT * FROM tasks WHERE task_id = %s"
            cursor.execute(sql, (task_id,))
            result = cursor.fetchone()
            if result:
                result = list(result)
                try:
                    if isinstance(result[3], str):
                        result[3] = json.loads(result[3])
                    if isinstance(result[4], str):
                        result[4] = json.loads(result[4])
                except json.JSONDecodeError:
                    pass  
                return tuple(result)
            return None
    finally:
        conn.close()

def register(user_name,password):
    conn = pymysql.connect(host=host, port=3306, user='sunpc', password='spc20020629', database='scseq')
    c = conn.cursor()

    c.execute('insert into users (name, password) values (%s, %s)',(user_name, password))
    conn.commit()

    c.close()
    conn.close()
    print("注册成功")

def login():
    conn = pymysql.connect(host=host, port=3306, user='sunpc', password='spc20020629', database='scseq')
    c = conn.cursor()
    user_name = 'spc'
    password = '123456'

    c.execute('select * from users where name = %s and password = %s', (user_name, password))
    result = c.fetchone()
    if result:
        print("登录成功")
    else:
        print("用户名不存在或密码错误")

# create_database()
# user_name = ''
# password = ''
# register(user_name,password)
