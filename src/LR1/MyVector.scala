package LR1;
/**
 * @author Ionkin Mikhail
 */
object MyVector{
  import math._
    def plus(a: Array[Double], b: Array[Double]): Array[Double] = {
      val c = new Array[Double](a.length); for (i <- 0 until a.length) c(i)=a(i)+b(i); c
    }
    def minus(a: Array[Double], b: Array[Double]): Array[Double] = {
      val c = new Array[Double](a.length); for (i <- 0 until a.length) c(i)=a(i)-b(i); c
    }
    def multi(a: Array[Double], f: Double):  Array[Double] = {
      Counter.countOpsMD+=a.length
      val c = new Array[Double](a.length); for (i <- 0 until a.length) c(i)=a(i)*f; c
    }
    def multi(a: Array[Double], b: Array[Double], i1: Int, i2: Int): Double = {
      Counter.countOpsMD+=(i2-i1)
      var sum = 0.0; for (i <- i1 until i2) sum+=a(i)*b(i); sum
    }
    def multi(a: Array[Double], B: Array[Array[Double]], j: Int, i1: Int, i2: Int): Double = {
      Counter.countOpsMD+=(i2-i1)
      var sum = 0.0; for (i <- i1 until i2) sum+=a(i)*B(i)(j); sum
    }
    def multi(a: Array[Double]): Double = {
      Counter.countOpsMD+=a.length
      var mult = 1.0; for (x <- a) mult*=x; mult
    }
    def indexOfMaxAbs(a: Array[Double], beginIndex: Int): Int = {
      var curInd = beginIndex; 
      for (i <- curInd+1 until a.length) if (abs(a(curInd)) < abs(a(i))) curInd = i
      curInd    
    }
    def multi(a: Array[Double], b: Array[Double]): Double = multi(a, b,0,a.length)
    def multi(a: Array[Double], B: Array[Array[Double]], j: Int): Double = multi(a,B,j,0,a.length)
    def normCubic(a: Array[Double]): Double = {var maxA = abs(a(0)); for (i <- 1 until a.length) maxA = max(maxA,abs(abs(a(i)))); maxA}
    def norm1(a: Array[Double]): Double = {
      var sum = 0.0; for(i <- 0 until a.length) sum+=abs(a(i)); sum
    }
    def norm2(a: Array[Double]) = sqrt(multi(a,a))
    def multiGM(a: Array[Double], b: Array[Double]): Array[Array[Double]] = {
      Counter.countOpsMD+=a.length*b.length
      val res = Array.ofDim[Double](a.length,b.length);
      for (i <- 0 until a.length) for (j <- 0 until b.length) res(i)(j) = a(i)*b(j)
      res
    }
    def toString(a: Array[Double], count: Int): String = {
      val format = "%."+count+'f'
      val out = new java.lang.StringBuilder()
      for (i <- 0 until a.length) out.append(format.format(a(i))).append(' ')
      out.toString()
    }
    def toString(a: Array[Double]): String = {//toString(a,2)
      val out = new java.lang.StringBuilder()
      for (i <- 0 until a.length-1) out.append("%.2f".format(a(i))).append(' ')
      out.append("%.2f".format(a(a.length-1)))
      out.toString()
    }
    def getMaxAndMin(a: Array[Double]): (Double,Double) = {
      var maxEV = abs(a(0)); var minEV = abs(a(0))
      for (i <- 0 until a.length){
        maxEV = max(maxEV,abs(a(i)))
        minEV = min(minEV,abs(a(i)))
      }
      (maxEV,minEV)
    }
    def indexOf[T: scala.reflect.ClassTag](a: Array[T], elem: T): Int = {
      for (i <- 0 until a.length) if (a(i).equals(elem)) return i
      -1
    }
    def swapInverse(xs: Array[Double],q:Array[Int]): Unit = {
      for (i <- 0 until xs.length) {
        val ind = indexOf(q,i)
        if (i!=ind) {
          val x = xs(i); xs(i)=xs(ind); xs(ind)=x;
          val y = q(i); q(i)=q(ind); q(ind)=y;
        }
      }        
    }
  }
