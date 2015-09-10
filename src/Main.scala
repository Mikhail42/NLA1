/**
 * @author Ionkin Mikhail
 */
object Main {
  import math._
  import LR1._
  /**
   * the relative error SLAE solutions using vector's norm "normCubic" (see object MyVector)
  */
  def error(a: Array[Double], b: Array[Double]): Double = {
    val c = MyVector.minus(a, b)
    MyVector.normCubic(c)/MyVector.normCubic(a)
  }
  def input() = MyMatrix.input(scala.io.Source.fromFile("/home/misha/Рабочий стол/matrix").mkString.split('\n'))
  def test(): Array[Array[Double]] = {
    import math._
    val n = (random*5.0).toInt+2
    val c = ((random+2.0)*100).toInt
    val A = Array.ofDim[Double](n,n)
    for (i <- 0 until n)
      for (j <- 0 until n) 
        A(i)(j) = c + Math.log1p(i*j)
    A
  }
  def main(args: Array[String]): Unit = {    
    def run(A: Array[Array[Double]]): Unit = {
      val n = A.length
      val x = new Array[Double](n)
      for (i <- 0 until n) x(i) = i+1
      val slae = new SLAE(A,x)
      //Устойчивость не гарантирована
      try{     
        Console.println(MyMatrix.toString(slae.LU_DWP.L))
        Console.println(MyVector.toString(slae.LU_DWP.xs))
        Console.println(MyMatrix.toString(slae.QR.Ak))
        Console.println(MyMatrix.toString(slae.QR.QR_HR._1))
        Console.println(MyMatrix.toString(slae.QR.QR_HR._2))
        Console.println(MyVector.toString(slae.QR.xs))
        Console.println(MyVector.toString(slae.QR.eigenvalues))
      }
      catch {
        case e: Exception => Console.println(e)
      }     
    }
    run(input())
    run(test())
  }
}
