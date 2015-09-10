package LR1;
/**
 * @author Ionkin Mikhail
 */
object MyMatrix{
  import math._
    def copy(inp: Array[Array[Double]], out: Array[Array[Double]]): Unit = {
      for (i <- 0 until inp.length) Array.copy(inp(i), 0, out(i), 0, inp.length)
    }
    def input(strs: Array[String]): Array[Array[Double]] = {
      val in1 = strs(0).split(' ').map{ x => x.toDouble }
      val n = in1.length
      val res =  Array.ofDim[Double](n,n);  
      Array.copy(in1,0,res(0),0,n);
      for (i <- 1 until n) res(i) = strs(i).split(' ').map{ x => x.toDouble }
      res
    }
    def toString(A: Array[Array[Double]]): String = {
      val out = new java.lang.StringBuilder()
      for (i <- 0 until A.length) out.append(MyVector.toString(A(i))).append('\n')
      out.toString()
    }
    def multi(A: Array[Array[Double]], f: Double): Array[Array[Double]] = {
      val n = A.length
      Counter.countOpsMD+=n*n
      val res = Array.ofDim[Double](n,n)
      for (i <- 0 until n) for (j <- 0 until n) res(i)(j) = A(i)(j)*f
      res
    }
    def multi(A: Array[Array[Double]], b: Array[Double]) : Array[Double] = {
      Counter.countOpsMD+=A.length*A.length
      val n = A.length
      val c = new Array[Double](n)
      for(i <- 0 until n) c(i) = MyVector.multi(A(i),b)
      c
    }
    def detDiagMatr(A: Array[Array[Double]]): Double = {
      Counter.countOpsMD+=A.length
      var res = 1.0; for (i <- 0 until A.length) res*=A(i)(i); res 
    }
    /**
     * p-norm, where p = infinity
     * the maximum absolute row sum of the matrix
     */
    def normMS(a: Array[Array[Double]]): Double = {
      import MyVector.norm1
      var maxA = norm1(a(0)); for (i <- 1 until a.length) maxA = max(maxA,norm1(a(i))); maxA
    }
    /**
     * p-norm, where p = 1
     * the maximum absolute column sum of the matrix
     */
    def normMC(a: Array[Array[Double]]): Double = {
      val n = a.length
      val mat = Array.ofDim[Double](n,n)
      for (i <- 0 until n) Array.copy(a, 0, mat, 0, n)
      transp(mat)
      normMS(mat)
    }
    /**
     * (p-norm, where p = 2)
     */
    def normSpectr(eigenvalues: Array[Double]): Double = {
      Counter.countOpsMD+=1
      val maxmin = MyVector.getMaxAndMin(eigenvalues)
      sqrt(maxmin._1/maxmin._2)
    }
    def multi(A: Array[Array[Double]], B: Array[Array[Double]]): Array[Array[Double]] = {
      val n = A.length
      val C = Array.ofDim[Double](n,n)
      for (i <- 0 until n)
        for (j <- 0 until n)
          C(i)(j) = MyVector.multi(A(i), B, j)
      C
    }
    def swapColumn(mat: Array[Array[Double]], i1: Int, i2: Int): Unit = {
      if (i1!=i2) 
        for (i <- 0 until mat.length){ 
          mat(i)(i1)+=mat(i)(i2); mat(i)(i2)=mat(i)(i1)-mat(i)(i2); mat(i)(i1)-=mat(i)(i2);
        }
    }
    def transp(mat: Array[Array[Double]]): Unit = {
      for (i <- 0 until mat.length)
        for (j<- 0 until i) {
          val x = mat(i)(j);
          mat(i)(j) = mat(j)(i)
          mat(j)(i) = x
        }
    }
  }
