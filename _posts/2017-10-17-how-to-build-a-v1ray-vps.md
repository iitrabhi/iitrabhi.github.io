---
layout: post
title: "如何搭建一个 V2Ray Server"
categories: tech
author: "Yizheng Huang"
meta: "Springfield"
---

打开 v2ray-core Github 的release界面，会看到开发者的一句话：“特别的爱给特别的你”，所以我沿用了这句话作为这篇服务器搭建教程的前言，也算是对开发者的致敬吧。
那么，什么是V2Ray呢？
由于天朝的封锁，科学上网变得越来越困难，V2Ray就是一个代理的平台，相比Shadowsocks它有更高的稳定性和加密的Vmess传输协议，被封锁的可能性大大低于Shadowsocks。具体关于它的介绍可以参考V2Ray的[官方文档和说明](https://www.v2ray.com/)
好的，闲话不说，直接进入正题

### 本文目录：
#### 1.新的服务器的一些必须的配置
- 准备
- 安装Git
- 配置终端fish

#### 2.如何在海外的服务器上搭建V2Ray
- 安装V2Ray
- 设置服务器时间
- 配置防火墙
- 设置服务端配置
- 设置服务器管理工具

####  3. 如何配置客户端
- Android 客户端的配置
- macOS 客户端的配置
- Windows 客户端的配置
- iOS 客户端的配置

---
### 1.新的服务器的一些必须的配置
#### 准备：
> - 一台位于海外的，不被封锁的服务器
> - 一个ssh软件，macOS 和 Linux 下一般使用ssh,Windows推荐使用 Xshell
> - 一个能够编辑json文件的文本编辑器(推荐Atom 和 SublimeText)

在这里我使用了通过Github Education申请到的[Github Student Pack](https://education.github.com/pack)，只要通过它的学生认证，你就可以得到一个开发工具包，其中里面有[DigitalOcean](https://www.digitalocean.com/)服务器的50美金代金券，DigitalOcean的服务器位于海外，不被封锁，非常适合开发者使用。具体怎么申请这个开发工具包和购买服务器，可以参考网上教程，我们这里不予赘述。

>实例教程配置环境
> - SanFrancisco 2区服务器
> - 最低价格5刀一个月的配置
> - Debian 9.2 x64


#### 安装Git：
Git是很重要的依赖，一般在第一次登录服务器的时候，就要配置好Git。
在登录DigitalOcean的服务器以后，我们首先要更新apt的源：
```bash
sudo apt-get update
```
更新完包管理器的源以后，直接下载安装Git:
```bash
sudo apt-get install git-core
```
这样我们就配置好了Git:
![配置Git](https://i.loli.net/2017/10/17/59e59dc46a152.png)
#### 配置终端fish：
fish 是一款很方便的终端shell 软件，个人喜好，如果不喜欢也可以直接就使用系统默认的bash
首先，下载安装fish shell：
```bash
sudo apt-get install fish
```
安装好了以后我们推荐fish的一个开源配置`Oh-My-Fish`，具体可以查看它的[官网](https://github.com/oh-my-fish/oh-my-fish)
打开终端，输入:
```bash
sudo apt-get install curl
```
安装好了curl以后：
```bash
curl -L https://get.oh-my.fish | fish
```
![安装好fish](https://i.loli.net/2017/10/17/59e59dff7df96.png)
这样就把fish完整地安装好了，接下来修改登录使用的shell：
输入以下命令来查看fish的安装位置:
```bash
which fish
```
然后再输入
```bash
chsh -s /你的fish的位置/
```
更改默认登录终端为fish
![更改默认shell](https://i.loli.net/2017/10/17/59e59e7fba8ce.png)
这样我们的开始前的设置就配置好了

### 2.如何在海外的服务器上搭建V2Ray
#### 安装V2Ray：
查看V2Ray用户手册，发现官网直接提供了下载的脚本，打开终端输入：
```bash
$ wget https://toutyrater.github.io/install-release.sh
```
执行一般不会出错，执行完成以后再次输入安装命令：
```bash
$ sudo  bash install-release.sh 
```
输出如下表示安装成功了:
![安装](https://i.loli.net/2017/10/17/59e5a06d8e277.png)
![安装成功](https://i.loli.net/2017/10/17/59e5a0942ffbd.png)
#### 设置服务器时间：
V2Ray对于服务器时间的设置要求十分严格，请求时间要是和服务器时间不一致（相差1min以内），就会认为是无效的服务器请求，所以安装好v2ray以后第一步就是要设置服务器实践与我们本地的时间一致，打开两个终端，一个连接到服务器另外一个保持本地。同时输入：
```bash
$ date -R
```
![观察时间](https://i.loli.net/2017/10/17/59e5a17f6aa38.png)
观察时间的差别，将服务器同步到本机时间即可，输入以下命令修改时间为东8区：
```bash
$ sudo  rm /etc/localtime
```
删除本地时间区文件，链接新的东八区时间文件，这里输入的是上海的：
```bash
$ sudo ln -s /usr/share/zoneinfo/Asia/Shanghai /etc/localtime
```
再次输入`date -R`检查时间：
![修改时区](https://i.loli.net/2017/10/17/59e5a294c67d3.png)
修改完时间以后，我们要同步硬件时间，这样重新启动服务器的时候就不会自动重置为原来的时间了，在这里我们要使用`ntpdate`这个命令，如果系统没有的话，自己下载一个：
```bash
$ sudo apt-get install ntpdate
```
接下来设置时间校准：
```bash
$ sudo ntpdate time.nist.gov
```
这样，我们关于时间的设置就结束啦
#### 配置防火墙：
这里我们使用的是德班系统，和ubuntu一样，用户可以选择在系统上设置防火墙：
```bash
$ sudo apt-get install ufw 
```
安装好ufw以后，我们将它启动，并且查看它的运行状态：
```bash
$ sudo ufw enable | status  
```
我们可以控制端口的开放和关闭，我们开放连接服务器用的22端口和后面要设置的代理端口，这个端口号可以自己设置：
```bash
$ ufw allow 22 
```
![防火墙](https://i.loli.net/2017/10/17/59e5a47f14c7e.png)
这样我们就配置好了防火墙
#### 设置服务端配置：
V2Ray的服务端设置是基于json的配置文件：`config.json`，一般放置在`/etc/v2ray/`目录下。这里不赘述关于它的各种配置，我们推荐新手直接使用一些服务端管理工具，如：[v2ray.fun](https://github.com/FunctionClub/v2ray.fun)
安装v2ray.fun命令：
```bash
$ wget -N --no-check-certificate https://raw.githubusercontent.com/FunctionClub/v2ray.fun/master/install.sh && bash install.sh
```
接下来我们稍微介绍一下关于服务端最基本的配置要求，安装好v2ray以后，在`/etc/v2ray/`目录下会有一个`config.json`文件,这个就是我们的代理服务端配置文件了，使用终端编辑软件打开，有这么一段位置：
```json
....省略一段....
 "inbound": {
        "port": 12345,
        "protocol": "vmess",
        "settings": {
            "clients": [
                {
                    "id": "5da8771d-1c5f-4007-a36a-1924e7cd25b5",
                    "level": 1,
                    "alterId": 64
                }
            ]
        },
        "streamSettings": {
            "network": "kcp"
        },
        "detour": {
            "to": "vmess-detour-565571"
        }
    },
....省略一段....
```
这里就是你的用户连接设置了，`port`代表的是连接的远程端口，就是我们在防火墙里面开启的那个端口，这个是可以由用户自己设置的，`protocol`代表着传输协议是`vmess`；`clients`里面包含的就是具体的用户连接设置；`id`在这里是用uuid来表示的唯一链接密码，可以去很多随机生成网站生成一个uuid，客户端的配置要和服务端一致才能够成功连接。同样`alterId`指的是一个动态的用户连接数，这个数越大越不容易给流量检测工具检测到，但是太大的话会严重占用计算机内存，所以一般我们设置一个在50到100之间的数字，注意这个参数在用户的配置文件里面也要保持一致。
同时，我们要设置v2ray的服务器随着服务器的开机自动启动，因为ssh关闭以后，服务器上的v2ray也就随之关闭了，我们只能通过一些开机配置脚本来保持它在后台的工作：
```bash
#!/bin/sh
### BEGIN INIT INFO
# Provides:          v2ray
# Required-Start:    $network $local_fs $remote_fs
# Required-Stop:     $remote_fs
# Default-Start:     2 3 4 5
# Default-Stop:      0 1 6
# Short-Description: socksv5 based proxy written by go.
# Description:       v2ray is a socksv5 proxy written by go. Connection can be crypto by aes or
#            des. this might help for people in China to corss GFW.
### END INIT INFO

# Author: Shell Xu <shell909090@gmail.com>
# Modify: Isulew Li <netcookies@gmail.com>

# PATH should only include /usr/* if it runs after the mountnfs.sh script
PATH=/sbin:/usr/sbin:/bin:/usr/bin  
DESC=v2ray             # Introduce a short description here  
NAME=v2ray             # Introduce the short server's name here  
DAEMON=/usr/bin/v2ray  #这里改成v2ray程序的绝对路径
PIDFILE=/var/run/$NAME.pid  
LOGFILE=/var/log/$NAME.log  
SCRIPTNAME=/etc/init.d/$NAME

DAEMON_OPTS="-config /etc/v2ray/config.json" #这里改成配置文件绝对路径

# Exit if the package is not installed
[ -x $DAEMON ] || exit 0

# Read configuration variable file if it is present
[ -r /etc/default/$NAME ] && . /etc/default/$NAME

# Load the VERBOSE setting and other rcS variables
. /lib/init/vars.sh

# Define LSB log_* functions.
# Depend on lsb-base (>= 3.0-6) to ensure that this file is present.
. /lib/lsb/init-functions

#
# Function that starts the daemon/service
#
do_start()  
{
    # Return
    #   0 if daemon has been started
    #   1 if daemon was already running
    #   2 if daemon could not be started
    #   3 if configuration file not ready for daemon
    start-stop-daemon --start --quiet --pidfile $PIDFILE --exec $DAEMON --test > /dev/null 
        || return 1
    start-stop-daemon --start --quiet --pidfile $PIDFILE --exec $DAEMON --background 
         --no-close -m -- $DAEMON_OPTS >> $LOGFILE 2>&1 
        || return 2
    chmod -f 600 $LOGFILE
    # Add code here, if necessary, that waits for the process to be ready
    # to handle requests from services started subsequently which depend
    # on this one.  As a last resort, sleep for some time.
}

#
# Function that stops the daemon/service
#
do_stop()  
{
    # Return
    #   0 if daemon has been stopped
    #   1 if daemon was already stopped
    #   2 if daemon could not be stopped
    #   other if a failure occurred
    start-stop-daemon --stop --quiet --retry=TERM/30/KILL/5 --pidfile $PIDFILE
    RETVAL="$?"
    [ "$RETVAL" = 2 ] && return 2
    # Wait for children to finish too if this is a daemon that forks
    # and if the daemon is only ever run from this initscript.
    # If the above conditions are not satisfied then add some other code
    # that waits for the process to drop all resources that could be
    # needed by services started subsequently.  A last resort is to
    # sleep for some time.
    start-stop-daemon --stop --quiet --oknodo --retry=0/30/KILL/5 --exec $DAEMON
    [ "$?" = 2 ] && return 2
    # Many daemons don't delete their pidfiles when they exit.
    rm -f $PIDFILE
    return "$RETVAL"
}

#
# Function that sends a SIGHUP to the daemon/service
#
do_reload() {  
    #
    # If the daemon can reload its configuration without
    # restarting (for example, when it is sent a SIGHUP),
    # then implement that here.
    #
    start-stop-daemon --stop --signal 1 --quiet --pidfile $PIDFILE
    return 0
}

case "$1" in  
  start)
    [ "$VERBOSE" != no ] && log_daemon_msg "Starting $DESC " "$NAME"
    do_start
    case "$?" in
        0|1) [ "$VERBOSE" != no ] && log_end_msg 0 ;;
        2) [ "$VERBOSE" != no ] && log_end_msg 1 ;;
    esac
  ;;
  stop)
    [ "$VERBOSE" != no ] && log_daemon_msg "Stopping $DESC" "$NAME"
    do_stop
    case "$?" in
        0|1) [ "$VERBOSE" != no ] && log_end_msg 0 ;;
        2) [ "$VERBOSE" != no ] && log_end_msg 1 ;;
    esac
    ;;
  status)
       status_of_proc "$DAEMON" "$NAME" && exit 0 || exit $?
       ;;
  reload|force-reload)
    #
    # If do_reload() is not implemented then leave this commented out
    # and leave 'force-reload' as an alias for 'restart'.
    #
    log_daemon_msg "Reloading $DESC" "$NAME"
    do_reload
    log_end_msg $?
    ;;
  restart|force-reload)
    #
    # If the "reload" option is implemented then remove the
    # 'force-reload' alias
    #
    log_daemon_msg "Restarting $DESC" "$NAME"
    do_stop
    case "$?" in
      0|1)
        do_start
        case "$?" in
            0) log_end_msg 0 ;;
            1) log_end_msg 1 ;; # Old process is still running
            *) log_end_msg 1 ;; # Failed to start
        esac
        ;;
      *)
        # Failed to stop
        log_end_msg 1
        ;;
    esac
    ;;
  *)
    #echo "Usage: $SCRIPTNAME {start|stop|restart|reload|force-reload}" >&2
    echo "Usage: $SCRIPTNAME {start|stop|status|reload|restart|force-reload}" >&2
    exit 3
    ;;
esac

```
复制这段配置，将它保存为v2ray，通过ssh上传到服务器上的`/etc/init.d/`目录下：
![上传自启动配置](https://i.loli.net/2017/10/17/59e5a8d6bbb5b.png)
为上传的文件加载权限：
```bash
$ sudo chmod +x /etc/init.d/v2ray
```
设置开机启动：

```bash
$ sudo update-rc.d v2ray defaults
```
这样就算配置好了
#### 设置服务器管理工具：
我们刚才安装好了`v2ray.fun`这个服务器管理脚本，安装完成以后，可以直接输入`v2ray`打开。
我们现将系统的v2ray服务开启：
```bash
$ sudo systemctl start v2ray
```
然后查看v2ray的运行状态：
```bash
$ service v2ray status
```
![v2ray状态](https://i.loli.net/2017/10/17/59e5aad0a26bf.png)
然后打开v2ray.fun管理工具：
![v2ray.fun](https://i.loli.net/2017/10/17/59e5ab3da50cd.png)
按照他的指示，可以很好的修改服务器配置：
![v2ray.fun修改配置](https://i.loli.net/2017/10/17/59e5ab0c21702.png)
到这里，我们的服务端配置就已经配置好了。

### 3. 如何配置客户端
配置好了我们的服务器端，现在很重要的就是客户端的配置，只有配置好了我们的客户端，我们才能科学上网。
#### Android 客户端的配置：

__准备__

- 客户端手机软件下载：[Actinium下载](https://github.com/V2Ray-Android/Actinium/releases/download/0.10.2/app-universal-release_aligned_signed.apk)
- Android 手机一台

__使用方法__

- 打开Actinium软件
- 点击软件右上角三个点，选择添加配置
- 下载上面配置，并将配置导入安卓手机
- 在软件中添加导入的配置文件，弹出是否转化的对话框，选择确定
- 点击右下角红色的按钮开始学习马克思主义
- 再次点击（打勾的形状）结束马克思主义的学习

#### macOS 客户端的配置：

### 准备
- 一台运行macOS的电脑
- 上述格式为：json的配置文件
- v2ray-core macOS版本[下载](https://github.com/v2ray/v2ray-core/releases/download/v2.41/v2ray-macos.zip)
- V2RayX 软件 [V2RayX.dmg下载](https://github.com/Cenmrev/V2RayX/releases/download/v0.7.8/V2RayX.dmg)
- HomeBrew macOS的包管理工具
- 终端软件
- 浏览器代理插件（这里使用Chrome 的SwitchyOmega插件）

### 使用方法
- 打开终端输入命令添加v2ray-core的homebrew Tap:
```bash
brew tap kofj/v2ray
```
- 使用brew 在系统安装v2ray-core:
```bash
brew install v2ray-core
```
__1. 使用图形界面客户端__

- 下载V2RayX软件并安装
- 下载上述格式为json的配置文件，并将文件名称设置为：`config.json` 
- 将`config.json`配置文件放置到目录：`/Applications/V2RayX.app/Contents/MacOS/`下
- 如果不想手动移动配置文件，可以采用图形化设置

__2. 使用V2Ray核心__

- 使用终端打开v2ray-core的话，直接将配置文件命名为：`config.json`放置到`/Users/用户名`下
- 然后直接在终端输入：
```bash
v2ray
```

__3. 上述两种方法都要使用的浏览器配置__

- 打开浏览器，保证之前开启的v2ray的终端窗口在回话不被关闭（或者是V2RayX软件不关闭）
- 在谷歌商店里面搜索：`SwitchyOmega` ,[下载地址](https://github.com/FelisCatus/SwitchyOmega/releases)并安装
- 打开浏览器`SwitchyOmega`设置：新建一个情景模式：
> 网址协议:默认
> 代理协议：SOCKS5
> 代理服务器：127.0.0.1 
> 代理端口：1080
> 删除“不代理的地址列表”里面原有的东西
- 设置好后，点击左下角：`应用选项`，并且将默认情景模式切换为自己新建的这个

__4. 能否正常使用检测__

- 打开一个新的终端：输入以下命令检测：
```bash
curl -v --socks5-hostname 127.0.0.1:1080 https://www.google.com/
```
- 若是显示大量谷歌官网数据，那么说明代理成功，去试一试浏览器的吧！

__注意__

- 使用v2ray-core不用打开V2RayX这个软件，但是得保持开启v2ray的终端一直运行
- 使用v2ray-core每次要使用代理的时候，只需要终端输入:
```bash
v2ray
```
#### Windows 客户端的配置：

__准备__
- 一台运行windows的计算机
- v2ray-core windows版本
- config.json 配置文件（本教程已贴出）

__使用方法__

__1. 安装v2ray-core__

- 安装v2ray-core，先下载官方软件[（下载x64）](https://github.com/v2ray/v2ray-core/releases/download/v2.41/v2ray-windows-64.zip)[（下载x32）](https://github.com/v2ray/v2ray-core/releases/download/v2.41/v2ray-windows-32.zip)
- 将下载好的文件解压，里面有一堆零散的文件。将自己的配置命名为`config.json`并且放在与解压得到的`v2ray`同一目录下，覆盖原有的`config.json`
- 将解压得到的文件夹整个添加到环境变量里面（具体怎么操作请自行搜索）

__2. 设置浏览器代理__

- 打开浏览器，保证之前开启的v2ray的终端窗口在回话不被关闭（或者是V2RayX软件不关闭）
- 在谷歌商店里面搜索：`SwitchyOmega` ,[下载地址](https://github.com/FelisCatus/SwitchyOmega/releases)并安装
- 打开浏览器`SwitchyOmega`设置：新建一个情景模式：
> 网址协议:默认
> 代理协议：SOCKS5
> 代理服务器：127.0.0.1 
> 代理端口：1080
> 删除“不代理的地址列表”里面原有的东西
- 设置好后，点击左下角：`应用选项`，并且将默认情景模式切换为自己新建的这个

__3. 设置完成__

- 设置完成以后每次要使用代理时，只要打开cmd或者按下`win+R`里面输入`v2ray`,正常运行v2ray即可
- 注意要保持运行v2ray的窗口一直在后台运行

#### iOS 客户端的配置：

别废话，先下载客户端`ShadowRocket`，现在这款软件已经在国区下架，有能力使用美区账号下载的童鞋或是能够找到下载账号的童鞋，我觉得配置已经不是什么难事了.....
