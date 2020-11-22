---
layout: post
title: "RxJava 引发的思考"
categories: tech
author: "Yizheng Huang"
meta: "Springfield"
---

在没有接触过 RxJava 以前，处理 Android 程序多线程是一件很伤脑的事情，曾经就有朋友向我抱怨：自己的Android 应用程序里面有大把大把的异步线程导致自己的代码可读性相当差劲，连自己都读不懂了，更何况日后自己去修改了。对此我也深有体会，由于自己对异步线程的管理不当，许多应用程序的线程出现了堵塞，死锁，甚至有的应用因为线程问题导致手机出现了 ANR 现象，简直了...

直到接触了 RxJava ，我才感受到原来异步线程的管理也可以这么简洁和有力，当然，这里所说的“简洁”并不是单纯的指代码的数量很少，而是在整体上有一个形象而清晰的思路，并且随着异步线程的增多和代码工程规模的增大，这种与传统异步写法和线程切换相比的优势更加突出。

我们来看一下以往传统异步线程的书写方式:

```java
new Thread() {
        @Override public void run () {
            super.run();
            for (File folder : folders) {
                File[] files = folder.listFiles();
                for (File file : files) {
                    if (file.getName().endsWith(".png")) {
                        final Bitmap bitmap = getBitmapFromFile(file);
                        getActivity().runOnUiThread(new Runnable() {
                            @Override
                            public void run() {
                                imageCollectorView.addImage(bitmap);
                            }
                        });
                    }
                }
            }
        }
    }.start();
```

这是一个向 imageCollectorView 添加本地图片的异步线程，现在我们再来看看这同样的代码，使用RxJava的写法：


```java
Observable.from(folders)
    .flatMap(new Func1<File, Observable<File>>() {
        @Override
        public Observable<File> call(File file) {
            return Observable.from(file.listFiles());
        }
    })
    .filter(new Func1<File, Boolean>() {
        @Override
        public Boolean call(File file) {
            return file.getName().endsWith(".png");
        }
    })
    .map(new Func1<File, Bitmap>() {
        @Override
        public Bitmap call(File file) {
            return getBitmapFromFile(file);
        }
    })
    .subscribeOn(Schedulers.io())
    .observeOn(AndroidSchedulers.mainThread())
    .subscribe(new Action1<Bitmap>() {
        @Override
        public void call(Bitmap bitmap) {
            imageCollectorView.addImage(bitmap);
        }
    });
```

有人会说：兄弟，你这个代码明明变多了啊！没错，就像之前说的一样，RxJava 侧重的是逻辑上面的简洁，你会发现在 RxJava 里面，每个方法所体现的思路都很清楚，整个异步操作的代码几乎就是一条直线下来的，一旦你的一步操作多了以后，RxJava 仍然可以一如既往的保持逻辑上的清晰。那么这是为什么呢？

其实，学习 RxJava，最最重要的就是理解 RxJava 里面蕴含的一种“观察者模式”，何为“观察者模式”呢？其实可以举一个例子：这个屋子里有一盏台灯，你走进这个屋子里，发现这个台灯时亮着的。

好了，在这个 case 里面，你就是观察者，而台灯就是动作的发生者。当然，这个故事还没有结束……你坐在屋子里静静地享受台灯的灯光（好吧，很少人会这么做），忽然，你发现台灯灭了（有可能是灯泡烧了），然后你采取了一系列措施让台灯又亮了起来。这个其实就是一个反馈的过程，台灯通过熄灭灯光告诉你：“我灭灯了”，然后你接收到了台灯给你发出的这个消息，并且采取了措施。

这是个很现实的例子，其实和 RxJava 的工作原理是一样的，在 RxJava 观察者模式里面，观察者，被观察者是最基本的组成成分。

现在，我们来看看怎么创建一个被观察者（Observable）：

```java
 Observable switcher=Observable.create(new Observable.OnSubscribe<String>(){

            @Override
            public void call(Subscriber<? super String> subscriber) {
                subscriber.onNext("On");
                subscriber.onNext("Off");
                subscriber.onNext("On");
                subscriber.onNext("On");
                subscriber.onCompleted();
            }
        });
        
```

当然，创建方式还不只这个,通过这种简单的just方式也可以创建Observable:

```java

Observable switcher=Observable.just("On","Off","On","On");

```
这样我们的一个观察者就创建好了，并且在被观察者发出的事件里面，我们定义了“On”、“Off”。
然后我们来继续创建被观察者(Subscriber)：

```java
 Subscriber light=new Subscriber<String>() {
            @Override
            public void onCompleted() {
                //被观察者的onCompleted()事件会走到这里;
                Log.d(TAG,"结束观察...\n");
            }

            @Override
            public void onError(Throwable e) {
                    //出现错误会调用这个方法
            }
            @Override
            public void onNext(String s) {
                //处理传过来的onNext事件
                Log.d(TAG,"handle this---"+s)
            }           
```
我们也可以这样稍微偷懒一下，通过一个人畜无害的Action：

```java
  Action1 light=new Action1<String>() {
                @Override
                public void call(String s) {
                    Log.d(TAG,"handle this---"+s)
                }
            }          
```

这样，我们的主角们就创建好啦！这个时候，有计算机基础的人就会问：观察者和被观察者在代码里面是通过怎样的一种方式建立联系的呢？其实吧，RxJava 的要素还有一个，那就是“订阅”。而相对应上述例子的“订阅”其实就是让观察者和被观察者建立一个联系的过程，就像你买杂志订阅一样，在代码的事件里面，观察者和被观察者之间的联系其实并没有系现实生活之中那么的复杂和多样，之间的联系全靠订阅来完成。

创建订阅的关系：

```java
switcher.subscribe(light);
```
这样，这个主角以及主角之间的相互关系都有了，我们再把代码整合一下：

```java
//创建被观察者，是事件传递的起点
Observable.just("On","Off","On","On")
        //这就是在传递过程中对事件进行过滤操作
         .filter(new Func1<String, Boolean>() {
                    @Override
                    public Boolean call(String s) {
                        return s！=null;
                    }
                })
        //实现订阅
        .subscribe(
                //创建观察者，作为事件传递的终点处理事件    
                  new Subscriber<String>() {
                        @Override
                        public void onCompleted() {
                            Log.d(TAG,"结束观察...\n");
                        }

                        @Override
                        public void onError(Throwable e) {
                            //出现错误会调用这个方法
                        }
                        @Override
                        public void onNext(String s) {
                            //处理事件
                            Log.d(TAG,"handle this---"+s)
                        }
        );
```

然而我们都知道，假如你在外面玩耍，然后有人忽然打电话告诉你你家着火了，需要马上回来或者打119救急（mdzz），那么那个打电话来告诉你这个消息的人其实在 RxJava 里面就是所谓的“订阅”成分，但是我们都知道，这样随随便便打莫名其妙的电话的人有可能是骗子，这个消息未必就是真的，如果说他有一定真实的成分，那也可能是添油加醋，经过处理加工的消息（其实是你家的台灯灭了而已）。

这个case在 RxJava 里面对应的就是：在订阅关系里面，可以对被观察者传出的事件进行各种处理和操作

在对应的编码操作里面，我们可以完成类似于Map变换还有FlatMap操作（and so on……） 
举个实例:

```java
//创建被观察者，获取所有班级
 Observable.from(getSchoolClasses())
                .flatMap(new Func1<SingleClass, Observable<Student>>() {
                    @Override
                    public Observable<Student> call(SingleClass singleClass) {
                        //将每个班级的所有学生作为一列表包装成一列Observable<Student>，将学生一个一个传递出去
                        return Observable.from(singleClass.getStudents());
                    }
                })
                .subscribe(
                //创建观察者，作为事件传递的终点处理事件    
                  new Subscriber<Student>() {
                        @Override
                        public void onCompleted() {
                            Log.d(TAG,"结束观察...\n");
                        }

                        @Override
                        public void onError(Throwable e) {
                            //出现错误会调用这个方法
                        }
                        @Override
                        public void onNext(Student student) {
                            //接受到每个学生类
                            Log.d(TAG,student.getName())
                        }
                    );
```

这样下来，你应该也发现了RxJava 之所以可以那么好的Manage所有的异步线程的原因了吧，代码的逻辑比传统的异步写作方式简洁得不知道哪里去了！

这个时候，RxJava 露出了一个狡猾的微笑：我厉害的东西你还没看到呢！没错，RxJava 厉害的不仅仅可以再订阅过程中处理一系列的异步操作和简洁的观察者模式，它还可以一句话切换任务发生的线程，实现线程的随时切换和高效的代码管理。

就像这样：

```java
          //new Observable.just()执行在新线程
  Observable.create(new Observable.just(getFilePath())
           //指定在新线程中创建被观察者
          .subscribeOn(Schedulers.newThread())
          //将接下来执行的线程环境指定为io线程
          .observeOn(Schedulers.io())
            //map就处在io线程
          .map(mMapOperater)
            //将后面执行的线程环境切换为主线程，
            //但是这一句依然执行在io线程
          .observeOn(AndroidSchedulers.mainThread())
          //指定线程无效，但这句代码本身执行在主线程
          .subscribeOn(Schedulers.io())
          //执行在主线程
          .subscribe(mSubscriber);
```
线程调度其实只有这两个方法：

> - subscribeOn（）它指示Observable在一个指定的调度器上创建（只作用于被观察者创建阶段）。只能指定一次，如果指定多次则以第一次为准
> - observeOn（）指定在事件传递（加工变换）和最终被处理（观察者）的发生在哪一个调度器。可指定多次，每次指定完都在下一步生效。

一句话，动作简洁有力。

这样，我们就把 RxJava 的背后的工作原理和基础用法理了一遍。记得当时我看完网络上的各种教程，内心是……完全蒙蔽的，因为在理解 RxJava 的过程之中，我感觉甚至要用到哲学的思维来看待这些代码。为什么它会比其他的代码更加简洁有力？为什么它会使用一种“观察者模式”作为写法的思路？随着对它的不断地了解，我感觉这个不仅仅是一些计算机代码，它的背后甚至可以拓展到其他学科的一种辩证的思维，就像著名的物理学家薛定谔关于猫的实验，打开盒子的一瞬间，决定猫死活的因素取决于你的观察，而你观察又是为了知道盒子里猫的死活，在 RxJava 里面就好像在订阅的这个关系里面,你给被观察者发出的信息一个决定性反馈一样。

这样的辩证思考不会被局限于技术领域，就拿Google最近开源的人工智能操作库 TensorFlow 来说，以张量流来完成机器的一个深度的学习，这个就涉及到了一个深度的反馈，代码在有了输出以后根据输出来重构代码，不断地实现学习和自我完善。拿“观察者模式”来解释，就是机器代码既是观察者，又是被观察者，而这个反馈递归的过程，就是日月自新的一个订阅关系。将代码的思维拓展到自然世界，根据环境的变化大自然自我调节，实现动态的平衡，并且具有一定的稳定性和抗力。这样的一个过程，其实已经说明了一点：计算机算法最后的归宿就是无限的贴近自然界，计算机语言无限地贴近自然语言，从而实现智能化。

啊，感觉自己扯了老远，RxJava 的使用其实远远不止这些，在处理观察事件的时候还可以有很多原本传统管理异步线程所没有的一些优点，这里就不说啦！时间不早了，赶紧收工回宿舍洗衣服，明天还有可怕的模电实验。。。

> 延伸阅读:
> - [关于RxJava最友好的文章](http://www.jianshu.com/p/6fd8640046f1)
> - [观察者模式的极致，RxJava运行流程及源码简析](http://www.jianshu.com/p/8091ade1e302)
> - [是时候学习RxJava了](http://www.jianshu.com/p/8cf84f719188)
> - [Rxjava学习：初识Rxjava](http://www.jianshu.com/p/264a1322669a)

