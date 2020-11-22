---
layout: post
title: "Bulid Java Project with Gradle"
categories: tech
author: "Yizheng Huang"
meta: "Springfield"
---


Gradle是Android Studio的默认依赖管理以及项目构建工具。所以，作为Java的初学者，应该怎么试着用Gradle来管理Java依赖以及构建Java项目呢？
学习Gradle 有助于将来Android等软件开发，深入理解Android一起系列编译器工作原理，这里只介绍利用Gradle构建基础的java项目，构建完成以后可以导入java IDE进行内部开发。

## Gradle的安装
Gradle 的安装要求系统本身具有__java__的开发环境
在安装Gradle之前可以先执行命令查看当前系统__java__版本：

```bash
$ java -version
java version "1.8.0_60"

```

### macOS (Homebrew) :

```bash
$ brew install gradle

```

### Linux（SDKMAN） :

```bash
$ sdk install gradle 3.4.1

```

### Windows (choco):


```bash
C:\> choco install gradle

```

## 使用 Gradle 构建基础项目

__Gradle__ 开发需要掌握一门新的脚本编程语言：__Groovy__
__Groovy__ 语言可以实现和 java 代码的相互调用，写法简单容易上手，这里不赘述，关于使用 __Gradle__ 构建基础 java 工程的方法如下：

### 新建目录，文件

手动建立一个目录，比如名叫：__gradle-java__ 然后进入该目录新建一个叫做 __build.gradle___的文件，这个是 __Gradle__ 的基本配置文件。


```bash
⋊> ~/Desktop mkdir gradle-java                                                                                                                                                         
⋊> ~/Desktop cd gradle-java/                                                                                                                                                           
⋊> ~/D/gradle-java nano build.gradle                                                                                                                                                   
⋊> ~/D/gradle-java

```
新建好build.gradle以后，打开build.gradle ，输入下面的关于构建java工程的配置信息：


```groovy
//bulid.gradle

apply plugin: 'java'

/*
 * 模板
 */
buildscript {
    repositories {
        maven {
            url 'http://dl.bintray.com/cjstehno/public'
        }
    }
    dependencies {
        classpath 'gradle-templates:gradle-templates:1.4.1'
    }
}
apply plugin:'templates'

```

上面的__buildscript__ 以及__apply plugin:'templates'__ 引入了一个叫做__templates的插件__，它能帮我们快速地创建项目所需要的文件。

然后，我们测试一下是否成功。

执行：

```bash
$ gradle -q tasks

```
发现控制台输出：

```bash
//省略一些输出……

Template tasks
--------------
createGradlePlugin - Creates a new Gradle Plugin project in a new directory named after your project.
createGroovyClass - Creates a new Groovy class in the current project.
createGroovyProject - Creates a new Gradle Groovy project in a new directory named after your project.
createJavaClass - Creates a new Java class in the current project.
createJavaProject - Creates a new Gradle Java project in a new directory named after your project.

//再省略一些输出……
```

这就说明成功了。

说明一下，__-q__参数是为了减少输出的某些无关紧要的东西，比如编译时间啊之类的，你省略了，也无关紧要。

### 从模板创建工程
完成了上面的工作后，当前的目录结构是这样子的。

```bash
gradlejava/
└── build.gradle

```

只有一个目录__gradlejava__,下面只有一个文件__build.gradle__
执行下面的命令：

```bash	
$ gradle -q initJavaProject

```

目录结构会变成下面这样子：

```bash
gradlejava/
├── LICENSE.txt
├── build.gradle
└── src
    ├── main
    │   ├── java
    │   └── resources
    └── test
        ├── java
        └── resources

```

然后执行：

```bash
$ gradle -q createJavaClass

```

按照提示输入类名,我输入的是 __cn.edu.sdu.online.hello.HelloWorld__ 其中__cn.edu.sdu.online__ 是山东大学学生在线的域名的倒写，__hello__ 是项目名称，__HelloWorld__ 是类名。我们再用同样的方法创建一个__cn.edu.sdu.online.hello.Greeter__

然后查看当下的目录结构：

```bash
gradlejava/
├── LICENSE.txt
├── build.gradle
└── src
    ├── main
    │   ├── java
    │   │   └── cn
    │   │       └── edu
    │   │           └── sdu
    │   │               └── online
    │   │                   └── hello
    │   │                       ├── Greeter.java
    │   │                       └── HelloWorld.java
    │   └── resources
    └── test
        ├── java
        └── resources

```

这就说明我们的java工程已经创建好了

### 写一个 HelloWorld 

我们开始在这个工程里面编写java代码：

```java
package cn.edu.sdu.online.hello;
 
public class HelloWorld {
    public static void main(String[] args) {
        Greeter greeter = new Greeter();
        System.out.println(greeter.sayHello());
    }
}

```
---

```java
package cn.edu.sdu.online.hello;

public class Greeter {
    public String sayHello() {
        return "安拉胡阿克巴";
    }
}

```

并且在__build.gradle__的结尾添加上:

```groovy
/*
 * Java 项目配置
 */
 
group = 'cn.edu.sdu.online'
apply plugin: 'application'
mainClassName = 'cn.edu.sdu.online.hello.HelloWorld'
repositories {
    mavenCentral()
}


```

（上面确定了程序的入口）
 然后我们执行下面的命令:

```bash
$ gradle -q run
 
 
```
 
 然后输出：
 
```bash
 安拉胡阿克巴
 
 
```
  
  证明HelloWorld写成功了，项目到此结束。
  
 ### 外部依赖
 
 在__bulid.gradle__末尾输入:
 
```groovy
 dependencies {
    compile "joda-time:joda-time:2.2"
}

  
```

我们就引入了__joda-time__这个库，需要时候，__gradle__会自动去下载的。不用我们操心。

然后我们修改__HelloWorld.java__

```java
package cn.edu.sdu.online.hello;
import org.joda.time.LocalTime;

 
public class HelloWorld {
    public static void main(String[] args) {
        LocalTime currentTime = new LocalTime();
        System.out.println("The current local time is: " + currentTime);
 
        Greeter greeter = new Greeter();
        System.out.println(greeter.sayHello());
    }
}

```

然后运行下面的命令:

```bash
$ gradle -q run


```

查看输出：

```bash
The current local time is: 18:14:19.243
安拉胡阿克巴


```

### 其他：
#### 程序打包 zip

```bash
$ gradle distzip


```

使用上述这个命令，你可以看到出现了一个 __build/distributions/gradlejava.zip__文件，我们就可以把这个文件给别人了。
他人收到zip文件，解压后，执行解压目录下的bin目录下的程序，就可以跑出结果了。

#### 程序打包 jar

要使用 gradle 打包 jar 包，得在 build.gradle 里面增加一条 task ，代码如下：

```gradle

jar {
    manifest {
        attributes "Main-Class": "$mainClassName"
    }

    from {
        configurations.compile.collect { it.isDirectory() ? it : zipTree(it) }
    }
}

dependencies {
    compile 'ch.qos.logback:logback-classic:1.1.2'
}


```
其中 __"jar"__ 是这个__task__的名称，在终端里面执行 gradle 编译打包的命令：

```bash
$ gradle jar

```

输出如下语句代表打包成功，如果build.gradle里面含有没有下载的依赖，gradle会自动下载，并且完成打包

```bash

mike@MikedeiMac ~/D/gradle-java> gradle jar
The Task.leftShift(Closure) method has been deprecated and is scheduled to be removed in Gradle 5.0. Please use Task.doLast(Action) instead.
        at build_blfxo9nnxbr4yoouhh2lfxoar.run(/Users/mike/Desktop/gradle-java/build.gradle:60)
Download https://repo1.maven.org/maven2/ch/qos/logback/logback-classic/1.1.2/logback-classic-1.1.2.pom
Download https://repo1.maven.org/maven2/ch/qos/logback/logback-parent/1.1.2/logback-parent-1.1.2.pom
Download https://repo1.maven.org/maven2/ch/qos/logback/logback-core/1.1.2/logback-core-1.1.2.pom
Download https://repo1.maven.org/maven2/org/slf4j/slf4j-api/1.7.6/slf4j-api-1.7.6.pom
Download https://repo1.maven.org/maven2/org/slf4j/slf4j-parent/1.7.6/slf4j-parent-1.7.6.pom
Download https://repo1.maven.org/maven2/ch/qos/logback/logback-classic/1.1.2/logback-classic-1.1.2.jar
Download https://repo1.maven.org/maven2/ch/qos/logback/logback-core/1.1.2/logback-core-1.1.2.jar
Download https://repo1.maven.org/maven2/org/slf4j/slf4j-api/1.7.6/slf4j-api-1.7.6.jar
:compileJava
:processResources NO-SOURCE
:classes
:jar

BUILD SUCCESSFUL

Total time: 9.088 secs

```

打包完成以后，就会在项目工程目录下面生成一个__build__文件夹，里面的__libs__就有编译打包好的 jar 包，在终端下面使用 jar 包，就可以直接运行刚才所编写的 java 程序了。

```bash

⋊> ~/D/g/b/libs java -jar gradle-java.jar                                                                                                                                              17:43:22
Hello World!

```

#### 清理编译结果

```bash
$ gradle clean


```

编译的文件就清理了。
这样，我们使用Gradle 来构建基础的java工程的步骤就讲完了。




