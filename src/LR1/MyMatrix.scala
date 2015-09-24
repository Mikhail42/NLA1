package LR1;
/**
 * @author Ionkin Mikhail
 */
object MyMatrix{
    import math._
    val THREAD_COUNT_SHIFT = 3
    val THREAD_COUNT = 1<<THREAD_COUNT_SHIFT
    def copy(inp: Array[Array[Double]], out: Array[Array[Double]], size: Int): Unit = {
      Counter.countIA+=size; Counter.countC+=(size+size)
      for (i <- 0 until size) Array.copy(inp(i), 0, out(i), 0, size)
    }
    def copy(inp: Array[Array[Double]], out: Array[Array[Double]]): Unit = copy(inp,out,inp.length)
    def input(strs: Array[String]): Array[Array[Double]] = {
      val in1 = strs(0).split(' ').map{ x => x.toDouble }
      val n = in1.length
      val res =  Array.ofDim[Double](n,n);  
      Array.copy(in1,0,res(0),0,n);
      for (i <- 1 until n) res(i) = strs(i).split(' ').map{ x => x.toDouble }
      Counter.countIA+=(n+1+1+n*n+n+n); Counter.countC+=n+n*n+n+n
      res
    }
    def toString(A: Array[Array[Double]]): String = {
      val out = new java.lang.StringBuilder()
      for (i <- 0 until A.length) out.append(MyVector.toString(A(i))).append('\n')
      Counter.countC+=A.length
      out.toString()
    }
    def multi(A: Array[Array[Double]], f: Double, size: Int): Array[Array[Double]] = {
      val n = size
      val res = Array.ofDim[Double](n,n)
      for (i <- 0 until n) for (j <- 0 until n) res(i)(j) = A(i)(j)*f
      Counter.countIA+=((n*n)<<1)+2; Counter.countMD+=n*n; Counter.countC+=((n*n)<<1)+2
      res
    }
    def multi(A: Array[Array[Double]], f: Double): Array[Array[Double]] = multi(A,f,A.length)
    def multi(A: Array[Array[Double]], b: Array[Double], size: Int) : Array[Double] = {
      val c = new Array[Double](size)
      for(i <- 0 until size) c(i) = MyVector.multi(A(i),b,0,size)
      Counter.countIA+=size*size+size; Counter.countC+=size
      c
    }
    def multi(A: Array[Array[Double]], b: Array[Double]) : Array[Double] = multi(A,b,A.length)
    def detDiagMatr(A: Array[Array[Double]]): Double = {
      Counter.countMD+=A.length; Counter.countC+=A.length; Counter.countIA+=A.length+1
      var res = 1.0; for (i <- 0 until A.length) res*=A(i)(i); res 
    }
    /**
     * p-norm, where p = infinity
     * the maximum absolute row sum of the matrix
     */
    def normMS(A: Array[Array[Double]], size: Int): Double = {
      import MyVector.norm1
      Counter.countC+=size+size; Counter.countIA+=size+1
      var maxA = norm1(A(0),size); for (i <- 1 until size) maxA = max(maxA,norm1(A(i),size)); maxA
    }
    def countZeros(A: Array[Array[Double]]): Int = {
      Counter.countC+=(A.length*A.length)<<1;  Counter.countIA+=A.length*A.length+1
      var count = 0; for (vect <- A) for (elem <- vect) if (MyDouble.isZero(elem)) count+=1; count;
    }
    def normMS(a: Array[Array[Double]]): Double = normMS(a,a.length);
    /**
     * p-norm, where p = 1
     * the maximum absolute column sum of the matrix
     */
    def normMC(A: Array[Array[Double]], size: Int): Double = {
      import MyVector.norm1
      Counter.countC+=(size-1)<<1; Counter.countIA+=size
      var maxA = norm1(A,0,size); for (i <- 1 until size) maxA = max(maxA,norm1(A,i,size)); maxA
    }
    def normMC(A: Array[Array[Double]]): Double = normMC(A,A.length);
    /**
     * (p-norm, where p = 2)
     */
    def normSpectr(eigenvalues: Array[Double]): Double = {
      Counter.countMD+=1; Counter.countIA+=1;
      val maxmin = MyVector.getMaxAndMin(eigenvalues)
      sqrt(maxmin._1/maxmin._2)
    }
    /**
     * size(A) = (n,n) = size(B), size in \doubleZ
     */
    def multi(A: Array[Array[Double]], B: Array[Array[Double]], size: Int): Array[Array[Double]] = {
      val n = size
      Counter.countC+=(n*n)<<1; Counter.countIA+=2+((n*n)<<1);
      val C = Array.ofDim[Double](n,n)
      for (k <- 0 to n/THREAD_COUNT)
        new Thread { 
          for (i <- THREAD_COUNT*k until min(THREAD_COUNT*(k+1),n))  
            for (j <- 0 until n)
              C(i)(j) = MyVector.multi(A(i), B, j, 0, n)
        }
      C
    }
    /**
     * size(A) = (n,n) = size(B)
     * multi(A,B,A.length)
     */
    def multi(A: Array[Array[Double]], B: Array[Array[Double]]): Array[Array[Double]] = MyMatrix.multi(A,B,A.length);
    def swapColumn(mat: Array[Array[Double]], i1: Int, i2: Int): Unit = {
      Counter.countC+=1; 
      if (i1!=i2) {
        Counter.countC+=mat.length; Counter.countIA+=mat.length*3; Counter.countPM+=mat.length*3
        for (i <- 0 until mat.length){ 
          mat(i)(i1)+=mat(i)(i2); mat(i)(i2)=mat(i)(i1)-mat(i)(i2); mat(i)(i1)-=mat(i)(i2);
        }
      }
    }
    def transp(mat: Array[Array[Double]]): Unit = {
      Counter.countC+=mat.length*mat.length; Counter.countIA+=3*mat.length*mat.length
      for (k <- 0 to mat.length/THREAD_COUNT)
        new Thread { 
          for (i <- THREAD_COUNT*k until min(THREAD_COUNT*(k+1),mat.length))  
            for (j<- 0 until i) {
              val x = mat(i)(j);
              mat(i)(j) = mat(j)(i)
              mat(j)(i) = x
            }
        }
    }
    def inverseTriangleUp(U: Array[Array[Double]]): Array[Array[Double]] = {
      val n = U.length
      Counter.countC += n*n+n+((n*(n+1))>>1); Counter.countMD += n+((n*(n+1))>>1); Counter.countIA += n*n+n+((n*(n+1))>>1)
      val res = Array.ofDim[Double](n,n)
      //[http://alexandr4784.narod.ru/kaplan5/kaplan_5_06.pdf p. 146] 
      for (i <- 0 until n) res(i)(i) = 1.0/U(i)(i);
      for (i <- 0 until n)  for (j <- i+1 until n) 
        res(i)(j) = -res(j)(j)*MyVector.multi(res(i), U, j, 0, j)
      res
    }
    def inverseTriangleDown(L: Array[Array[Double]]): Array[Array[Double]] = {
      val n = L.length
      Counter.countC += n*n+n+((n*(n+1))>>1); Counter.countMD += n+((n*(n+1))>>1); Counter.countIA += n*n+n+((n*(n+1))>>1)
      val res = Array.ofDim[Double](n,n)
      //[http://alexandr4784.narod.ru/kaplan5/kaplan_5_06.pdf p. 149] 
      for (i <- 0 until n) res(i)(i) = 1.0/L(i)(i);
      for (i <- 0 until n) for (j <- 0 until i) 
        res(i)(j) = -res(i)(i)*MyVector.multi(L(i),res, j, 0, i)
      res
    }
  }
