/**
 * @author Ionkin Mikhail
 */
object Main {
  import math._
  import LR1._
  def input() = MyMatrix.input(scala.io.Source.fromFile("/home/misha/Рабочий стол/matrix").mkString.split('\n'))
  def test(): Array[Array[Double]] = {
    import math._
    val n = (random*10.0).toInt+3
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
      def xsJoS() = slae.METHOD_ID match{
        case 1 => slae.JacobiAndSeidel.xsJacobi
        case 2 => slae.JacobiAndSeidel.xsSeidel
      }
      Console.println("Исходная матрица:")
      Console.println(MyMatrix.toString(A))
      //Устойчивость не гарантирована
      try{     
        val out1 = new java.lang.StringBuilder(
            " Таблица 1. \n 1. Метод ").append(
            slae.METHOD_ID match{case 1 => "Якоби"; case 2 => "Зейдаля";}).append(
            "\n 2. Относительная погрешность решения: "+slae.error(x,xsJoS())).append(
            "\n 3. Число операций: ").append(
            Counter.countOpsMD)
        
        Console.println(out1)
        Console.println()
      }
      catch {
        case e: Exception => Console.println(e)
      }  
      try{      
        val out2 = new java.lang.StringBuilder(
            " Таблица 2. \n 1. Определитель матрицы системы: ").append(
            slae.LU_DWP.det).append(
            "\n 2. Макс. и мин. собственные числа матрицы A (не A*A^T): ").append(
            {val mm = MyVector.getMaxAndMin(slae.QR.eigenvalues); (mm._1.toString()+" и "+mm._2.toString())}).append(
            "\n 3. Число обусловленности матрицы системы: ").append(
            slae.LU.cond)
        Console.println(out2)
        Counter.countOpsMD = 0
        Console.println()
        Console.println()
      } catch {
        case e: Exception => Console.println(e)
      }  
    }
    run(input())
    run(test())
  }
}
