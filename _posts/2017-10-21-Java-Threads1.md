---
layout: post
title: "Java 并发编程手记 (1)"
categories: tech
author: "Yizheng Huang"
meta: "Springfield"
---

## 本文目录
- 故事缘由
- synchronized 关键字学习
  - 什么是 synchronized ?
  - synchronized 关键字的作用域
  - 为什么要使用 synchronized 关键字？
  - synchronized 关键字的特性
  - synchronized 同步锁使用的优化

## 故事缘由

前一阵面试完了今日头条的实习，技术面试的时候由于我在时间内提早写完了代码，所以面试官又加问了几个问题。其中一个是关于 Java 线程安全的，以前在工程上很少遇到这类问题，但是这类问题在面试的时候出现的概率是非常高的，题目不难，大概是这个样子：

```java
// 这段代码是否是线程安全的？如果不是，怎么修改？

public class A {
      int a = 0;

      void funA(){
           a++;
      }

      int funB(){
           int b = a + 1;
           return b;
      }
}
```
这个题目是很典型的并发问题，面试官叫我现场给他修改看看，于是我把方法改成了:

```java
public class A {
      int a = 0;

      synchronized void funA(){
           a++;
      }

      synchronized int funB(){
           int b = a + 1;
           return b;
      }
}
```
后来思考，感觉自己当时的回答其实不好，这两个方法里面，容易造成并发问题的其实只有一个变量 `a` ，考虑到程序的效率， `synchronized` 关键字实际上是很重的，对于这种情况可以直接使用 Java 里面自带的一些原子类，在这里可以使用 `AtomicInteger`，改进的代码如下：

```java
public class A {
      AtomicInteger a = new AtomicInteger(0);

      void funA(){
           count.incrementAndGet();
      }

      int funB() {
           int b = a.addAndGet(1);
           return b;
     }
}

```

在面完这个问题之后，我觉得自己在 Java 并发编程上面的基础仍然薄弱，还需要系统的学习和提升，毕竟这个还是属于编程语言的基础，对于面试来说，即使其他的算法解答和项目介绍得非常好，侥幸通过了面试，但是基础不牢固还是没法在今后的工作和研究中走远。

这几日无聊，打算重新系统地对 Java 并发编程进行一个复习，同时记录在此也算当作个人的笔记，如果有不足或者是错误的地方，欢迎大家指正。

## synchronized 关键字学习
### 什么是 synchronized ?
在介绍这个关键字的时候，我想先说一个生活中的场景，假如有一个卫生间，一次只能一个人使用，如果两个人同时挤到卫生间里面就会出问题 (误)。有些代码也是一样，在多线程的环境下，如果多个线程同时调用了某段代码，这些代码处理的结果在不同线程里面就有可能出现不同步的问题。所以，我们使用了 `synchronized` ，它是一个同步锁，它的作用就是保证同一时间代码调用的同步性。

换句话来说，`synchronized` 就像是卫生间包厢前面的锁，一个人进去了以后，他拿到了这把锁，把自己锁在包厢里面，这样同一个时间就只有他能够享用卫生间了，没有锁的其他人是无法访问这个卫生间的，就只能在卫生间门口排队，等待里面的那个人上完卫生间出来，释放锁，把锁交给下一个排队的人。

**那么，说了半天，什么是锁呢？（感觉很抽象）**
我们先来阅读一段代码：
```java
public class SynchronizedUse {
    private int count = 10;
    // 锁对象
    private Object o = new Object();

    public void m() {
        synchronized (o) { // 想要执行下面一段代码，必须先拿到 o 的锁
            count--;
            System.out.println(Thread.currentThread().getName() + " count= " + count);
        }
    }
}
```
在这一段代码里面，我们新建了一个 Object 对象，并且使用 `synchronized` 作用在了了这个对象 `o` 上。很多人理解锁概念的时候出现了误解，认为 `synchronized` “锁”住的是 `synchronized` 作用的代码块，其实这是不对的，`synchronized` 锁着的是对象。想要执行上述代码中被 `synchronized` 修饰的代码块，只有拿到对象`o`的锁才行。

由于对象只有一个，所以代码保证了每次只会有一个线程能够拿到这把锁，只有一个线程能够执行被锁着的代码。

### synchronized 关键字的作用域
> 那我是不是每次想要同步一段代码，都得新建一个 Object 对象呢？

不是的，阅读下面这段代码会发现，其实我们可以直接使用 `this` 对象来进行代码的锁定，这是一种简化的写法。

```java
public class SynchronizedUse02 {
    private int count = 10;

    private void m() {
        synchronized (this) {
            count--;
            System.out.println(Thread.currentThread().getName() + " count= " + count);
        }
    }
}
```
`synchronized` 关键字除了使用花括号`{}`修饰一个代码块（同步代码块）以外，还可以直接修饰一个方法，被修饰的方法也被称为同步方法，其实和直接修饰一段代码块相同，修饰方法锁定的也是 `this` 对象，如下面代码所示，这个写法等同于`SynchronizedUse02 `(上一个demo)。 

```java
public class SynchronizedUse03 {

    private int count = 10;

    public synchronized void m() { // 等同于 synchronized(this){ ...}
        count--;
        System.out.println(Thread.currentThread().getName() + " count= " + count);
    }
}
```

同时，`synchronized` 关键字还可以修饰一个静态的方法，其作用的范围是整个静态方法，作用的对象是这个类的所有对象。了解 Java `static` 关键字的同学们都知道，作用在 `static` 方法上面的 `synchronized` 关键字实际上没有作用在任何实例化的对象上，而是直接作用在类对象上面。如下面代码所示：

```java
public class SynchronizedUse04 {

    private static int count = 10;

    // synchronized 使用在静态方法上的时候，相当于锁定了 class
    public synchronized static void m() {
        count--;
        System.out.println(Thread.currentThread().getName() + " count= " + count);
    }

    // 相当于这个方法
    public static void mm() {
        // 实际上是反射
        synchronized (SynchronizedUse04.class) {
            count--;
        }
    }
}
```
这么理解，即使我实例化了不同的`SynchronizedUse04`对象，在不同的线程里面调用静态方法 `m()` ，它仍然会保持同步，因为静态方法是属于类的，而不是属于对象的，它对该类的所有对象都会保持同步。

### 为什么要使用 synchronized 关键字？

之前说了那么多，我们一直都在强调一个"同步"，那么为什么在高并发程序中，同步那么重要呢，我们用一个小的demo来说明:

```java
public class SynchronizedUse05 implements Runnable {

    private int count = 10;

    // 如果不加锁，那么容易出现重复的数字,且得不到顺序打印的数字。
    // 每个 synchronized 的代码块，都代表一个原子操作，是最小的一部分，不可分。
    public synchronized void run() {
        count--;
        System.out.println(Thread.currentThread().getName() + " count = " + count);
    }

    /*
     * main 方法里面，实际上是一个对象只启动了一个方法。
     * 但是在 for 循环里面新建了多个线程来访问一个对象 t 。
     * */
    public static void main(String[] args) {
        SynchronizedUse05 t = new SynchronizedUse05();
        for (int i = 0; i < 8; i++) {
            // 新建的8个线程都去访问 t 里面的 run() 方法。
            new Thread(t, "THREAD" + i).start();
        }
    }
}
```

我们在这个例子里面启动了8个不同的线程，每个线程都会调用`t`里面的`run()`方法，而每个线程都对变量`count`进行了自减操作。源代码为这个类的`run()`方法加了同步锁，如果把锁去掉，我们就很容易在每次之中得到不同的运行结果，或者说出现重复的数字，且每次运行得不到顺序打印的数字，如下图所示，我们可以这样理解：
![不加同步锁的后果](https://i.loli.net/2019/09/07/zaM5e3AoYwZWHmU.png)

`Thread 1` 和 `Thread 2` 几乎同时执行代码，它们都拿到了值为10的`count`并对其进行了修改，所以这两个线程就会输出一样的`count`，之后由于线程对于资源的抢占式得到，所以陆陆续续输出结果的线程也不会是按照顺序的，这也是为什么`Thread 1`执行完毕以后`Thread 4`接下去执行的原因了，而加了锁之后可以消除这样的问题。


不加锁的运行结果 (每次未必一致)：

```
THREAD5 count = 4
THREAD6 count = 3
THREAD7 count = 2
THREAD3 count = 6
THREAD2 count = 7
THREAD0 count = 9
THREAD1 count = 8
THREAD4 count = 5
```

加了锁的运行结果 (运行结果每次一致):
```
THREAD0 count = 9
THREAD7 count = 8
THREAD6 count = 7
THREAD5 count = 6
THREAD4 count = 5
THREAD3 count = 4
THREAD2 count = 3
THREAD1 count = 2
```

同时，对多线程读写代码加上同步锁，还可以避免常见的 **“脏读问题” (Dirty Read)**，所谓脏读问题，指的是对业务写方法加锁，对业务读方法不加锁，那么同一个时间写入的线程只能有一个，但是读取却不受限制，这样，在写入的同时另外一个线程进行读取，就容易读取到错误的数据，轻则报出空指针异常，重则读取到匪夷所思的错误数据。我们可以通过下面的这个demo来体会一下：

```java
public class SynchronizedUse07 {
    String name;
    double balance;

    private synchronized void set(String name, double balance) {
        this.name = name;

        // 写入的时候加锁，要是写入时间中还在执行一些其他的程序，这时候读程序在另外一个线程中
        // 读取信息，写入工作还没完成，就容易读取到错误的信息。
        try {
            Thread.sleep(2000);
        } catch (InterruptedException e) {
            e.printStackTrace();
        }

        this.balance = balance;
    }

    private /* synchronized */ double getBalance(String name) {
        return this.balance;
    }

    public static void main(String[] args) {
        SynchronizedUse07 a = new SynchronizedUse07();

        new Thread(() -> a.set("zhangsan", 100.0)).start();

        try {
            TimeUnit.SECONDS.sleep(1);
        } catch (InterruptedException e) {
            e.printStackTrace();
        }

        System.out.println(a.getBalance("zhangsan"));

        try {
            TimeUnit.SECONDS.sleep(2);
        } catch (InterruptedException e) {
            e.printStackTrace();
        }

        System.out.println(a.getBalance("zhangsan"));
    }
}
```

在上述代码中，我们把读取账户余额的方法 ` getBalance(String name)` 上面的同步锁 `synchronized` 注释掉了，保留了写入（初始化）方法的同步锁，这样运行就会产生脏读问题。为了让代码问题突出，我们在`set()`里面加入了一段延时程序:

```java
     try {
            Thread.sleep(2000);
        } catch (InterruptedException e) {
            e.printStackTrace();
        }
```
在业务逻辑中，写入数据可能是一个耗时操作（这里我们使用了2秒的延时）,而如果读取操作不上锁，在写入操作还未完成的时候就开始读取，就会读取到错误数据（这里我们在延时1秒后开始读取，这时候写入操作尚未完成），所以第一次我们读取到的数据是 0.0（未初始化值），而当写入操作完成以后，我们再次进行读取操作，这时候才读取到正确的数据 100.0。

程序输出结果如下：

```java
0.0 // 1秒的时候开始读取，写入未完成
100.0 // 3秒的时候读取，写入（耗时约2秒）已完成
```

通过上面两个例子，我们现在对 `synchronized` 关键字的作用和用法有了一个初步的了解。

### synchronized 关键字的特性

> synchronized 关键字修饰的同步方法能不能和非同步方法同时调用？

从之前脏读问题的demo可以很轻易地得到答案：__可以__。
为此我们编写了一个demo，感兴趣的朋友可以尝试运行体会一下:

```java
public class SynchronizedUse06 {
    // m1 是同步方法，请问在执行m1 的过程之中，m2能不能被执行？ 回答：当然可以
    private synchronized void m1() {
        System.out.println(Thread.currentThread().getName() + " m1 start...");

        try {
            Thread.sleep(10000);
        } catch (InterruptedException e) {
            e.printStackTrace();
        }

        System.out.println(Thread.currentThread().getName() + " m1 end.");
    }

    private void m2() {
        try {
            Thread.sleep(5000);
        } catch (InterruptedException e) {
            e.printStackTrace();
        }

        System.out.println(Thread.currentThread().getName() + " m2.");
    }

    public static void main(String[] args) {
        SynchronizedUse06 t = new SynchronizedUse06();

//        new Thread(() -> t.m1(), "t1").start();
//        new Thread(() -> t.m2(), "t1").start();

        new Thread(t::m1, "t1:").start();
        new Thread(t::m2, "t2:").start();
    }
}
```
这里的`m1()`是一个同步方法，其中我们为它加入了一个长达10秒的延时，`m2()`是一个非同步方法。我们新建了两个线程，线程`t1`在运行方法`m1()`的时候，我们启动线程`t2`运行方法`m2()`，可以看到，`m1()`和`m2()`是可以同时运行的：

```java
t1: m1 start...
t2: m2.
t1: m1 end.
```
在启动线程这里我们使用了 Java8 的新特新：Lambda表达式，这样的作用主要是简化写法（语法糖），感兴趣的同学可以另外了解。

> 对于同步方法来说，它可以和非同步方法同时调用，那么对于同步方法之间的相互调用来说，这是否可行？

我们先来看一段小程序：
```java
/**
 * synchronized 关键字的使用08
 * 一个同步方法可以调用另外一个同步方法
 * 一个线程已经拥有某个对象的锁，再次申请的时候仍然会得到该对象的锁。
 * 也就是说 synchronized 获得的锁是可以重入的。
 *
 * @author huangyz0918
 */
public class SynchronizedUse08 {
    private synchronized void m1() {
        System.out.println("m1 start...");

        try {
            TimeUnit.SECONDS.sleep(1);
        } catch (InterruptedException e) {
            e.printStackTrace();
        }

        m2();
    }

    private synchronized void m2() {
        try {
            TimeUnit.SECONDS.sleep(1);
        } catch (InterruptedException e) {
            e.printStackTrace();
        }

        System.out.println("m2.");
    }

    public static void main(String[] args) {
        SynchronizedUse08 t = new SynchronizedUse08();
        t.m1();
    }
}
```

运行结果:

```java
m1 start...
m2.
```

在代码中，我们调用了同步方法`m1()`，`m2()`也得到了调用，说明同步方法直接是可以相互调用的。方法`m1()`和`m2()`需要的是同一把锁，所以当`m2()`在`m1()`中被调用的时候，`m1()`在持有锁的前提下再次申请获得这把锁来执行`m2()`，这是可以的，因为`synchronized`获得的锁是可以重入的。

> 那么对于继承的同步方法来说，子类方法是否也是同步方法？

__不是的，synchronized 关键字不能被继承。__

虽然可以使用`synchronized` 来定义方法，但 `synchronized` 并不属于方法定义的一部分，因此，`synchronized` 关键字不能被继承。如果在父类中的某个方法使用了 `synchronized` 关键字，而在子类中覆盖了这个方法，在子类中的这个方法默认情况下并不是同步的，而必须显式地在子类的这个方法中加上`synchronized`关键字才可以。当然，还可以在子类方法中调用父类中相应的方法，这样虽然子类中的方法不是同步的，但子类调用了父类的同步方法，因此，子类的方法也就相当于同步了。具体可以见下面的 demo：

```java
public class SynchronizedUse09 {

    public static void main(String[] args) {
        T t = new T();
        t.m();
    }

    synchronized void m() {
        System.out.println("m start...");

        try {
            TimeUnit.SECONDS.sleep(1);
        } catch (InterruptedException e) {
            e.printStackTrace();
        }

        System.out.println("m end.");
    }
}

class T extends SynchronizedUse09 {
    @Override
    synchronized void m() {
        System.out.println("child m start...");
        super.m();
        System.out.println("child m end.");
    }
}
```
这里的`m()`显式地加上了关键字`synchronized`，所以它也是一个同步方法。如果不加这个关键字，`m()`就不会得到同步，当然，在调用到父类`super.m()`的时候，仍然需要获得锁才能够执行，否则方法将一直在`super.m()`上面等待。

那么对于一个程序，在获得了锁开始执行的时候，忽然在同步代码块里面出现了异常，跳出了同步代码，那么锁是否会释放？

> 程序在执行的过程中，如果出现异常，默认情况下锁会被释放。

我们可以试着建立一个 demo 研究一下：

```java
public class SynchronizedUse10 {
    int count = 0;

    synchronized void m() {
        System.out.println(Thread.currentThread().getName() + " start...");
        while (true) {
            count++;
            System.out.println(Thread.currentThread().getName() + " count = " + count);

            try {
                TimeUnit.SECONDS.sleep(1);
            } catch (InterruptedException e) {
                e.printStackTrace();
            }

            if (count == 5) {
                int i = 1 / 0; // 此处抛出异常，锁将被释放，若是不想锁被释放可以进行 catch 使循环继续。
            }
        }
    }

    public static void main(String[] args) {
        SynchronizedUse10 t = new SynchronizedUse10();
        new Thread(t::m, "t1").start();

        try {
            TimeUnit.SECONDS.sleep(3);
        } catch (InterruptedException e) {
            e.printStackTrace();
        }

        new Thread(t::m, "t2").start();
    }
}
```

运行结果：

```java
t1 start...
t1 count = 1
t1 count = 2
t1 count = 3
t1 count = 4
t1 count = 5
t2 start...
Exception in thread "t1" java.lang.ArithmeticException: / by zero
t2 count = 6
	at learn.multithreading.synchronizedlearn.SynchronizedUse10.m(SynchronizedUse10.java:32)
	at java.base/java.lang.Thread.run(Thread.java:844)
t2 count = 7
t2 count = 8
t2 count = 9
t2 count = 10
t2 count = 11
t2 count = 12
```

这里我们使用了一个死循环，在`count == 5`地时候，我们手动的抛出了一个异常，使得线程`t1`从同步方法`m()`中跳出，并且释放了锁，我们可以看到，在`t1`抛出异常以后，`t2`迅速获得了锁并且开始执行。说明程序在执行的过程中，如果出现异常，默认情况下锁会被释放。 所以，在并发处理的过程中，如果出现了异常一定要多加小心，不然可能会发生不一致的情况。比如在一个 Web App 处理的过程中，多个 servlet 线程共同访问一个资源，这时候如果异常处理不合适， 在第一个线程里面抛出异常，其他线程就会进入同步代码区，有可能会访问到异常时产生的数据。

### synchronized 同步锁使用的优化

虽然现在 Java 已经对 `synchronized` 锁进行了很大幅度的优化，不过相对与其他同步机制来说，`synchronized` 的效率还是比较低的，所以在使用这个关键字的时候，我们要注意锁的粒度，避免不必要的计算资源浪费。

举个例子：

```java
public class SynchronizedUse11 {
    private int count = 0;

    synchronized void m1() {
        try {
            TimeUnit.SECONDS.sleep(2);
        } catch (InterruptedException e) {
            e.printStackTrace();
        }

        count++; // 业务逻辑中只有这句话需要同步，这时候不需要给整个方法上锁。

        try {
            TimeUnit.SECONDS.sleep(2);
        } catch (InterruptedException e) {
            e.printStackTrace();
        }
    }

    void m2() {
        try {
            TimeUnit.SECONDS.sleep(2);
        } catch (InterruptedException e) {
            e.printStackTrace();
        }

        // 业务逻辑中只有这一段语句需要同步，这时不应该给整个方法上锁
        // 采用细粒度的锁，可以使线程争用的时间变短，从而提高效率。
        synchronized (this) {
            count++;
        }

        try {
            TimeUnit.SECONDS.sleep(2);
        } catch (InterruptedException e) {
            e.printStackTrace();
        }
    }
}
```

在这里，我们只有`count++` 不符合原子性，所以容易出现同步问题，没有必要对整个方法进行上锁，所以说`synchronized`关键字的使用还是很灵活的，在编码的时候要时时刻刻考虑到效率的问题。

同时，我们应该理解`synchronized`锁定的本质，`synchronized`其实锁定的是堆内存中的对象，所以当一个锁定的对象的属性发生了改变，或者说，锁定对象的引用指向了堆内存中的一个新的对象的时候，锁也会改变。在实际使用当中，我们应该避免发生这样的情况：

```java
public class SynchronizedUse12 {

    public static void main(String[] args) {
        SynchronizedUse12 t = new SynchronizedUse12();
        new Thread(t::m, "t1").start();

        try {
            TimeUnit.SECONDS.sleep(3);
        } catch (InterruptedException e) {
            e.printStackTrace();
        }

        Thread t2 = new Thread(t::m, "t2");

        t.o = new Object(); // 锁的对象改变，所以 t2 得以执行，不然 t2 永远得到不了执行的机会。
        t2.start();
    }

    Object o = new Object();

    void m() {
        synchronized (o) { // 本质上说明了: 锁的位置是锁在堆内存的对象上，而不是栈内存对象的引用里面。
            while (true) {
                try {
                    TimeUnit.SECONDS.sleep(2);
                } catch (InterruptedException e) {
                    e.printStackTrace();
                }

                System.out.println(Thread.currentThread().getName());
            }
        }
    }
}
```

输出结果：

```java
t1
t1
t2 // 从这里开始，t2线程也开始了运行
t1
t2
t1
t2
t1
....
```

在这个小程序里面，我们可以看到我们使用` t.o = new Object();`将锁的对象指向了堆内存里面的一个新的对象，这个时候线程`t2`其实已经和`t1`所需要的不是同一把锁了，所以线程`t2`也开始运行。不然处于死循环方法`m()`里面，不等到`t1`执行完成，`t2`永远也得不到执行的机会。

**最后，我们在实际的开发过程中，经常要尽力避免使用字符串对象进行锁定。**为什么呢？如果你使用了某个类库，它锁定了字符串对象 A ，这时候你在自己的源代码里面又锁定了字符串对象 A ，两段代码不经意之间使用了同一把锁，这样就会出现很诡异的死锁阻塞，而且对于你来说这样的问题很难被排查出来。并且对于字符串来说，两个相同的字符串其实指向的是同一个内存地址，所以看似使用的不是同一把锁，实际上不然：

```java
public class SynchronizedUse13 {

    private String s1 = "Hello";
    private String s2 = "Hello";

    // 两个字符串s1和s2实际上指向了同一个堆内存的对象
    // 栈内存里面存放原始变量和对象的引用句柄
    // 堆内存里面存放的是对象的实例

    void m1() {
        synchronized (s1) {
        }
    }

    void m2() {
        synchronized (s2) {
        }
    }
}
```

总之，Java 同步锁 `synchronized` 的使用是很灵活的，需要在实践中不断总结和反复记忆。

