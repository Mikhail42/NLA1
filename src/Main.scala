/**
 * @author Ionkin Mikhail
 */
object Main {
  import math._
  import LR1._
  def input() = MyMatrix.input(scala.io.Source.fromFile("/home/misha/Рабочий стол/matr").mkString.split('\n'))
  def testSimple(): Array[Array[Double]] = {
    import math._
    val n = 200
    val A = Array.ofDim[Double](n,n)
    for (i <- 0 until n)
      for (j <- 0 until n) 
        A(i)(j) = Math.sqrt(random*1000)
    A
  }
  def test(): Array[Array[Double]] = {
    import math._
    val n = (random*100.0).toInt+3
    val c = (random+2.0)*100.0
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
      val mySLAE = new SLAE(A,x)
      import slae._
      //Console.println("Исходная матрица:")
      //Console.println(MyMatrix.toString(A))
      Console.println(" Таблица 1.");
      Console.println("Метод \t\t\t\t Относительная погрешность решения\t Количество операций (mod 2^64)")
      try{     
        Console.print("1. LU-разложение\t\t ");
        mySLAE.setLU;  
        Console.print("%.18f".format(slae.ISLAE.error(x, mySLAE.getLU.xs))+"\t\t\t")
        Console.println(Counter)
        Counter.setZeros()
      } catch {
        case e: Exception => Console.println(e.getMessage)
      }
      try{
        Console.print("2. LU-разложение c выбором в. э. ");
        mySLAE.setLUP; 
        Console.print("%.18f".format(slae.ISLAE.error(x, mySLAE.getLUP.xs))+"\t\t\t")
        Console.println(Counter)
        Counter.setZeros()
      } catch {
        case e: Exception => Console.println(e.getMessage)
      }
      try{  
        Console.print("3. "+"Якоби"+"\t\t\t ");
        mySLAE.setJacobiAndSeidel; 
        Console.print({val s = "%.18f".format(slae.ISLAE.error(x,mySLAE.getJacobiAndSeidel.xsJoS())); (s+{if (s.length()<10) "\t\t" else ""})}+"\t\t\t")
        Console.println(Counter)
        Counter.setZeros()
      } catch {
        case e: Exception => Console.println(e.getMessage)
      }
      Console.println()
      Console.println(" Таблица 2.");
      try{      
        Console.print("1. Определитель матрицы системы: \t\t ") 
        Console.println(mySLAE.getLUP.det)
      } catch {
        case e: Exception => Console.println(e.getMessage)
      }
      try{
        Console.print("2. Макс. и мин. собственные числа матрицы A*A^T: ")
        Console.println({mySLAE.setQR; 
            val m1 = MyVector.getMaxAndMin(mySLAE.getQR.eigenvalues); (m1._1.toString()+" и "+m1._2.toString())
            //val m2 = MyVector.getMaxAndMin(mySLAE.getQR.eigenvalues); (m2._1.toString()+" и "+m2._2.toString())
            })
      } catch {
        case e: Exception => Console.println(e.getMessage)
      }
      try{
        Console.print("3. Число обусловленности матрицы системы: \t ")
        Console.println(mySLAE.getLU.getAndSetCond())
      } catch {
        case e: Exception => Console.println(e.getMessage)
      }
      Console.println("Число операций (mod 2^64, после последнего измерения): "+Counter.toString())
      Counter.setZeros()
      Console.println()
      Console.println()
    }
    run(testSimple())
    run(test())
  }
}
