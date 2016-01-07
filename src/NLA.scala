import structLA._
import slae.SLAE
import math._

object NLA {
  
 
  /**
   * @param args:
   *  args(0): String --fileNameWithParameters. default: /home/misha/Рабочий стол/param
   * fileNameWithParameters: 
   *  0: String -- inFileName                      (exp.: /home/misha/Рабочий стол/matr)
   *  1: toInt  -- QR_ITER                         (exp.: 20)
   *  2: toDouble -- QR_ERROR                      (exp.: 1e-3)
   *  3: toInt -- n (size of matrixs, n*n)         (exp.: 100)
   *  4: toInt isPrintMatrix                       (1 - yes, 0 - not)
   *  5: toInt isInput (=> run(input(inFileName))  (1 - yes, 0 - not)
   *  6: toInt isSimpleTest (=> run(testSimple(n)) (1 - yes, 0 - not)
   *  7: toInt isTest (=> run(test(n))             (1 - yes, 0 - not)
   *  8: normID                                  p of p-norms, where p = NORM_ID. If (p==3) => (p = Infinity)
   *  9: methodID                                Jacobi or Seidel 1 - Jacobi method 2 - Seidel method
   *  10: DOUBLE_ERROR                           (exp: 1e-10. exp., structLA.MyDouble.isZero(1.e-11)==true)
   */
  def main(args: Array[String]): Unit = { 

    slae.Function.METHOD_ID = 1     // methodID
    slae.Function.NORM_ID   = 1     // normID
    slae.Function.QR_ERROR  = 1e-3  // QR_ERROR
    slae.Function.QR_ITER   = 40    // QR_ITER
    MyDouble.doubleError = 1e-10 // DOUBLE_ERROR)
       
    val n = 200
    val A = testSimple(n)
    run(A)
    val B = test(n)
    run(B)
  }
    
  def run(A: Array[Array[Double]]): Unit = {    
      import slae._
      
      val n = A.length
      val x = new Array[Double](n)
      for (i <- 0 until n) x(i) = i+1
      val mySLAE = new SLAE(A,x)
      
      Console.println(" Таблица 1.");
      Console.println("Метод \t\t\t Относительная погрешность решения\t Количество операций (mod 2^64)")
      try{     
        Console.print("1. LU-разложение\t ")
        Console.print(slae.Function.error(x, mySLAE.lu.xs)+"\t\t\t")
        //Console.println(MyMatrix.toString(slae.ISLAE.getLLT(A)))
        //Console.println(MyMatrix.toString(slae.ISLAE.getLU(A, false, null)._1))
        //Console.println(MyMatrix.toString(slae.ISLAE.getLU(A, false, null)._2))
        Console.println(Counter)
        Counter.setZeros
      } catch {
        case e: Exception => Console.println(e.getMessage)
      }
      try{
        Console.print("2. LUP-разложение\t ")
        Console.print(slae.Function.error(x, mySLAE.lup.xs) + "\t\t\t")
        Console.println(Counter)
        Counter.setZeros
      } catch {
        case e: Exception => Console.println(e.getMessage)
      }
      try{  
        Console.print("3. "+"Якоби"+"\t\t ");
        Console.print(slae.Function.error(x, mySLAE.jacobiAndSeidel.xsJoS) +"\t\t\t")
        Console.println(Counter)
        Counter.setZeros
      } catch {
        case e: Exception => Console.println(e.getMessage)
      }
      Console.println()
      Console.println(" Таблица 2.");
      try{      
        Console.print("1. Определитель матрицы системы: \t\t ") 
        Console.println(mySLAE.lup.det)
      } catch {
        case e: Exception => Console.println(e.getMessage)
      }
      try{
        Console.print("2. Макс. и мин. собственные числа матрицы A*A^T: ")
        Console.println({val m1 = MyVector.getMaxAndMin(mySLAE.qr.eigenvalues); (m1._1.toString()+" и "+m1._2.toString())})
      } catch {
        case e: Exception => Console.println(e.getMessage)
      }
      try{
        Console.print("3. Число обусловленности матрицы системы: \t ")
        Console.println(mySLAE.lup.cond)
      } catch {
        case e: Exception => Console.println(e.getMessage)
      }
      Console.println("Число операций (mod 2^64, после последнего измерения): "+Counter.toString())
      Counter.setZeros
      Console.println()
      Console.println()
    }
  
  def input(name: String) = MyMatrix.input(scala.io.Source.fromFile(name).mkString.split('\n'))
  
  def testSimple(n: Int): Array[Array[Double]] = {
    import math._
    val A = Array.ofDim[Double](n,n)
    for (i <- 0 until n; j <- 0 until n) 
      A(i)(j) = Math.sqrt(random*1000)
    A
  }
  
  def test(n: Int): Array[Array[Double]] = {
    import math._
    val c = (random+2.0)*100.0
    val A = Array.ofDim[Double](n,n)
    for (i <- 0 until n; j <- 0 until n) 
      A(i)(j) = 1.0/(i+j+1.0)    // c + Math.log1p(i*j) 
    A
  }
}