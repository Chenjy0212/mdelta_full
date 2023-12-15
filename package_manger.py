import os

# 引用python模块库
'''
字典={库名: 安装情况}
0未安装，1已安装
'''

dict_pkg_name = {    # 数据处理库
                'ipywidgets': '0',  # mdetal_menu 显示窗口
                'munkres': '0',
                'matplotlib': '0',
                'networkx': '0',
                '--upgrade pip': '0'  # 更新包 # pip install --upgrade package_name
                }


# 国内镜像源
dict_pkg_NET = {'清华': R'pip3 install {pkg_name} -i https://pypi.tuna.tsinghua.edu.cn/simple ',
                '阿里云': R'pip3 install {pkg_name} -i http://mirrors.aliyun.com/pypi/simple/',
                '中国科技大学': R'pip3 install {pkg_name} -i https://pypi.mirrors.ustc.edu.cn/simple/',
                '华中理工大学': R'pip3 install {pkg_name} -i http://pypi.hustunique.com/',
                '山东理工大学': R'pip3 install {pkg_name} -i http://pypi.sdutlinux.org/',
                '豆瓣': R'pip3 install {pkg_name} -i http://pypi.douban.com/simple',
                '豆瓣2': R'pip install {pkg_name} -i http://pypi.douban.com/simple/ --trusted-host pypi.douban.com',
                '官网2': R'pip3 install {pkg_name}',
                '官网': R'pip install {pkg_name}'
                }

count = 3  # 设置循环运行3次
while count:
    try:

        print('import packages: Begin')

        import ipywidgets
        import munkres
        import matplotlib
        import networkx

        print('import packages: END')
        break

    except Exception as err:
        print('Import model error: ' + str(err))
        print('\n')
        # 遍历库名
        for key_pkg_name in dict_pkg_name:
            if dict_pkg_name[key_pkg_name] == '0':
                # 遍历网址进行安装
                for key_pkg_NET in dict_pkg_NET:
                    try:
                        pip_install = dict_pkg_NET[key_pkg_NET].format(pkg_name=key_pkg_name)
                        print('Installing: '+pip_install)
                        os.system(pip_install)
                        dict_pkg_name[key_pkg_name] = '1'
                        print('\n')
                        break
                    except Exception as err:
                        print('Fail ...: ' + str(err))
                        # print('\n')
                        pass
    finally:
        count -= 1


