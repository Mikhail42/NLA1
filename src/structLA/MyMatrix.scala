package structLA;

import scala.math.max
import scala.math.sqrt

/**
 * @author Ionkin Mikhail
 */
object MyMatrix{
    import math._
    
    def input(strs: Array[String]): Array[Array[Double]] = {
      val in1 = strs(0).split(' ').map{ x => x.toDouble }
      val n = in1.length
      val res =  Array.ofDim[Double](n,n);  
      Array.copy(in1,0,res(0),0,n);
      for (i <- 1 until n) res(i) = strs(i).split(' ').map{ x => x.toDouble }
      Counter.countIA+=(n+1+1+n*n+n+n); Counter.countC+=n+n*n+n+n
      res
    }
    
    def multi(A: Array[Array[Double]], f: Double, size: Int): Array[Array[Double]] = {
      val n = size
      val res = Array.ofDim[Double](n,n)
      for (i <- 0 until n) for (j <- 0 until n) res(i)(j) = A(i)(j)*f
      Counter.countIA+=((n*n)<<1)+2; Counter.countMD+=n*n; Counter.countC+=((n*n)<<1)+2
      res
    }
    
    def multi(A: Array[Array[Double]], f: Double): Array[Array[Double]] = multi(A,f,A.length)
    
    def copy(inp: Array[Array[Double]], out: Array[Array[Double]], size: Int): Unit = {
      Counter.countIA+=size; Counter.countC+=(size+size)
      for (i <- 0 until size) Array.copy(inp(i), 0, out(i), 0, size)
    }
    
    def copy(inp: Array[Array[Double]], out: Array[Array[Double]]): Unit = copy(inp,out,inp.length)
    
    def toString(A: Array[Array[Double]]): String = {
      val out = new java.lang.StringBuilder()
      for (i <- 0 until A.length) out.append(MyVector.toString(A(i))).append('\n')
      out.toString()
    }

    def multi(A: Array[Array[Double]], b: Array[Double], size: Int) : Array[Double] = {
      val c = new Array[Double](size)
      for(i <- 0 until size) c(i) = MyVector.multi(A(i),b,0,size)
      c
    }
    
    def multi(A: Array[Array[Double]], b: Array[Double]) : Array[Double] = multi(A,b,A.length)
    
    def detDiagMatr(A: Array[Array[Double]]): Double = {
      var res = 1.0; for (i <- 0 until A.length) res*=A(i)(i); res 
    }
    
    /**
     * p-norm, where p = infinity
     * the maximum absolute row sum of the matrix
     */
    def normMS(A: Array[Array[Double]], size: Int): Double = {
      import MyVector.norm1
      var maxA = norm1(A(0),size); 
      for (i <- 1 until size) maxA = max(maxA,norm1(A(i),size)); maxA
    }
    
    def countZeros(A: Array[Array[Double]]): Int = {
      var count = 0; for (vect <- A) for (elem <- vect) if (MyDouble.isZero(elem)) count+=1; count;
    }
    
    def normMS(a: Array[Array[Double]]): Double = normMS(a,a.length);
    
    /**
     * p-norm, where p = 1
     * the maximum absolute column sum of the matrix
     */
    def normMC(A: Array[Array[Double]], size: Int): Double = {
      import MyVector.norm1
      var maxA = norm1(A,0,size); for (i <- 1 until size) maxA = max(maxA,norm1(A,i,size)); maxA
    }
    
    def normMC(A: Array[Array[Double]]): Double = normMC(A,A.length);
    
    /**
     * (p-norm, where p = 2)
     */
    def normSpectr(eigenvalues: Array[Double]): Double = {
      val maxmin = MyVector.getMaxAndMin(eigenvalues)
      sqrt(maxmin._1/maxmin._2)
    }
    
    /**
     * size(A) = (n,n) = size(B), size in \doubleZ
     */
    def multi(A: Array[Array[Double]], B: Array[Array[Double]], size: Int): Array[Array[Double]] = {
      val n = size
      val C =  Array.ofDim[Double](n,n)
      for (i <- 0 until n) for (j <- 0 until n) C(i)(j) = MyVector.multi(A(i), B, j, 0, n)
      C
    }   
    /**
     * size(A) = (n,n) = size(B)
     * multi(A,B,A.length)
     */
    def multi(A: Array[Array[Double]], B: Array[Array[Double]]): Array[Array[Double]] = MyMatrix.multi(A,B,A.length);
    def swapColumn(mat: Array[Array[Double]], i1: Int, i2: Int): Unit = {
      if (i1!=i2) {
        for (i <- 0 until mat.length){ 
          mat(i)(i1)+=mat(i)(i2); mat(i)(i2)=mat(i)(i1)-mat(i)(i2); mat(i)(i1)-=mat(i)(i2);
        }
      }
    }
    def transp(mat: Array[Array[Double]]): Unit = {
      val n = mat.length
      for (i <- 0 until n)
            for (j<- 0 until i) {
              val x = mat(i)(j);
              mat(i)(j) = mat(j)(i)
              mat(j)(i) = x
            }
    }
    def inverseTriangleUp(U: Array[Array[Double]]): Array[Array[Double]] = {
      val n = U.length
      val res = Array.ofDim[Double](n,n)
      //[http://alexandr4784.narod.ru/kaplan5/kaplan_5_06.pdf p. 146] 
      for (i <- 0 until n) res(i)(i) = 1.0/U(i)(i);
      for (i <- 0 until n)  for (j <- i+1 until n) 
        res(i)(j) = -res(j)(j)*MyVector.multi(res(i), U, j, 0, j)
      res
    }
    def inverseTriangleDown(L: Array[Array[Double]]): Array[Array[Double]] = {
      val n = L.length
      val res = Array.ofDim[Double](n,n)
      //[http://alexandr4784.narod.ru/kaplan5/kaplan_5_06.pdf p. 149] 
      for (i <- 0 until n) res(i)(i) = 1.0/L(i)(i);
      for (i <- 0 until n) for (j <- 0 until i) 
        res(i)(j) = -res(i)(i)*MyVector.multi(L(i),res, j, 0, i)
      res
    }
  }