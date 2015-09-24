package LR1;
/**
 * @author Ionkin Mikhail
 * elementwise operations {.*,.+,.-,./} mean the same thing as in {Matlab,Scilab}
 */
object MyVector{
  import math._
    /** a(:).+b(:) */
    def plus(a: Array[Double], b: Array[Double]): Array[Double] = {
      Counter.countIA+=a.length+1; Counter.countC+=a.length; Counter.countPM+=a.length;
      val c = new Array[Double](a.length); for (i <- 0 until a.length) c(i)=a(i)+b(i); c
    }
    /** a(:).-b(:) */
    def minus(a: Array[Double], b: Array[Double]): Array[Double] = {
      Counter.countIA+=a.length+1; Counter.countC+=a.length; Counter.countPM+=a.length;
      val c = new Array[Double](a.length); for (i <- 0 until a.length) c(i)=a(i)-b(i); c
    }
    /** f*a(:) */
    def multi(a: Array[Double], f: Double):  Array[Double] = {
      Counter.countMD+=a.length; Counter.countIA+=a.length+1; Counter.countC+=a.length; 
      val c = new Array[Double](a.length); for (i <- 0 until a.length) c(i)=a(i)*f; c
    }
    /** a(i1:i2-1).*b(i1:i2-1) = sum(ai*bi : i in [i1,i2)) */
    def multi(a: Array[Double], b: Array[Double], i1: Int, i2: Int): Double = {
      Counter.countMD+=(i2-i1); Counter.countIA+=(i2-i1)+1; Counter.countC+=i2-i1; Counter.countPM+=i2-i1;
      var sum = 0.0; for (i <- i1 until i2) sum+=a(i)*b(i); sum
    }
    /** a(i1:i2-1).*B(i1:i2-1)(j) */
    def multi(a: Array[Double], B: Array[Array[Double]], j: Int, i1: Int, i2: Int): Double = {
      Counter.countMD+=(i2-i1); Counter.countIA+=(i2-i1)+1; Counter.countC+=i2-i1; Counter.countPM+=i2-i1;
      var sum = 0.0; for (i <- i1 until i2) sum+=a(i)*B(i)(j); sum
    }
    /** a(:).*a(:) */
    def multi(a: Array[Double]): Double = {
      Counter.countMD+=a.length; Counter.countIA+=a.length+1; Counter.countC+=a.length; 
      var mult = 1.0; for (x <- a) mult*=x; mult
    }
    /** 
     *  Search index maximal element in the abs(a(:)). Search begins with the beginIndex. For example: 
     *  @beginIndex = 1
     *  @param a (-11, -10, 9, 2)
     *  @return 1
     *  */ 
    def indexOfMaxAbs(a: Array[Double], beginIndex: Int): Int = {
      var curInd = beginIndex;
      Counter.countIA+=a.length-curInd; Counter.countC+=(a.length-curInd-1)<<1; 
      for (i <- curInd+1 until a.length) if (abs(a(curInd)) < abs(a(i))) curInd = i; curInd    
    }
    /** 
     *  inner product
     *  multi(a, b,0,a.length) 
     **/ 
    def multi(a: Array[Double], b: Array[Double]): Double = multi(a, b,0,a.length)
    /** multi(a,B,j,0,a.length) */
    def multi(a: Array[Double], B: Array[Array[Double]], j: Int): Double = multi(a,B,j,0,a.length)
    /** max(abs(a(0:n-1))) */
    def normCubic(a: Array[Double], n: Int) = {
      var maxA = abs(a(0)); for (i <- 1 until n) maxA = max(maxA,abs(a(i))); maxA
    }
    /** normCubic(a,a.length) */
    def normCubic(a: Array[Double]): Double = normCubic(a,a.length);
    /** sum(abs(a(0:n-1))) */
    def norm1(a: Array[Double], n: Int): Double ={var sum = 0.0; for(i <- 0 until n) sum+=abs(a(i)); sum}
    /** norm1(a,a.length) */
    def norm1(a: Array[Double]): Double = norm1(a,a.length)
    /** sum(abs(A(0:n-1)(j))) */
    def norm1(A: Array[Array[Double]],j: Int, n: Int): Double ={var sum = 0.0; for(i <- 0 until n) sum+=abs(A(i)(j)); sum}
    /** norm1(A,j,A.length) */
    def norm1(A: Array[Array[Double]],j: Int): Double = norm1(A,j,A.length)
    /** sqrt(multi(a,a,0,n)) */
    def norm2(a: Array[Double], n: Int) = sqrt(multi(a,a,0,n))
    /** norm2(a,a.length) */
    def norm2(a: Array[Double]): Double = norm2(a,a.length)
    /**
     * outer product
     *  a*b^T 
     **/
    def multiGM(a: Array[Double], b: Array[Double]): Array[Array[Double]] = {
      Counter.countMD+=a.length*b.length; Counter.countC+=a.length*b.length; Counter.countIA+=a.length*b.length+1;
      val res = Array.ofDim[Double](a.length,b.length);
      for (i <- 0 until a.length) for (j <- 0 until b.length) res(i)(j) = a(i)*b(j)
      res
    }
    /** 
     *  @count the number of digits after the decimal point. 
     *  format = '%'+@count.toString+'f'.  
     *  @return sum(format.format(a(:))+' ')
     *  */
    def toString(a: Array[Double], count: Int): String = {
      val format = "%."+count+'f'
      val out = new java.lang.StringBuilder()
      for (i <- 0 until a.length) out.append(format.format(a(i))).append(' ')
      out.toString()
    }
    /** toString(a,4) */
    def toString(a: Array[Double]): String = {//toString(a,4)
      val out = new java.lang.StringBuilder()
      for (i <- 0 until a.length-1) out.append("%.4f".format(a(i))).append(' ')
      out.append("%.4f".format(a(a.length-1)))
      out.toString()
    }
    /** 
     *  @return (max(abs(a)),min(abs(a)))
     *  */
    def getMaxAndMin(a: Array[Double]): (Double,Double) = {
      Counter.countC+=a.length; Counter.countIA+=(a.length<<1);
      var maxEV = abs(a(0)); var minEV = abs(a(0))
      for (i <- 1 until a.length){
        maxEV = max(maxEV,abs(a(i)))
        minEV = min(minEV,abs(a(i)))
      }
      (maxEV,minEV)
    }
    /**
     * It returns the index of the first matched element
     * @param a = (2, 3, 3, 4)
     * @param elem = 3
     * @return 1 
     */
    def indexOf[T: scala.reflect.ClassTag](a: Array[T], elem: T): Int = {
      Counter.countC+=a.length; 
      for (i <- 0 until a.length) if (a(i).equals(elem)) return i
      -1
    }
    /** rearranges a vector xs in accordance with the permutation vector q 
     *  @param x = (1,2)
     *  @param q = (1,0)
     *  @return (2,1)
     *  */
    def swapInverse(xs: Array[Double],q:Array[Int]): Unit = {
      Counter.countC+=(xs.length<<1); Counter.countIA+=xs.length*7;
      for (i <- 0 until xs.length) {
        val ind = indexOf(q,i)
        if (i!=ind) {
          val x = xs(i); xs(i)=xs(ind); xs(ind)=x;
          val y = q(i); q(i)=q(ind); q(ind)=y;
        }
      }        
    }
  }
