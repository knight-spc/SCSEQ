import socket

ipv4s = socket.gethostbyname_ex(socket.gethostname())[2]
print(ipv4s[0])
print("************************")


import platform
import sys

print("Python版本:", sys.version)
print("操作系统:", platform.system(), platform.release())
print("Python路径:", sys.executable)