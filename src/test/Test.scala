package test
/**
 * @author Ionkin Mikhail
 */
object Test {
  import math._
  import LR1._
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
   */
  def main(args: Array[String]): Unit = { 
    val fileNameWithParameters = if (!args.isEmpty) args(0) else "/home/misha/Рабочий стол/param"
    val params = scala.io.Source.fromFile(fileNameWithParameters).mkString.split('\n')
    
    val inFileName = params(0)
    val QR_ITER = params(1).toInt
    val QR_ERROR = params(2).toDouble
    val n = params(3).toInt
    val isPrintMatrix = params(4).toInt
    val isInput = params(5).toInt
    val isSimpleTest = params(6).toInt
    val isTest = params(7).toInt
    val normID = params(8).toInt
    val methodID = params(9).toInt
    
    if (isInput == 1) run(input(inFileName))
    if (isSimpleTest == 1) run(testSimple(n))
    if (isTest == 1) run(test(n))
    def run(A: Array[Array[Double]]): Unit = {
      val n = A.length
      val x = new Array[Double](n)
      for (i <- 0 until n) x(i) = i+1
      val mySLAE = new SLAE(A,x)
      slae.ISLAE.METHOD_ID = methodID
      slae.ISLAE.NORM_ID = normID
      import slae._
      if (isPrintMatrix == 1) {
        Console.println("Исходная матрица:")
        Console.println(MyMatrix.toString(A))
      }
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
        Console.print({val s = "%.18f".format(slae.ISLAE.error(x,mySLAE.getJacobiAndSeidel.xsJoS())); 
                    (s+{if (s.length()<10) "\t\t" else ""})}+"\t\t\t")
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
        mySLAE.setQR; mySLAE.getQR.ERROR = QR_ERROR;  mySLAE.getQR.ITER = QR_ITER
        Console.println({val m1 = MyVector.getMaxAndMin(mySLAE.getQR.eigenvalues); (m1._1.toString()+" и "+m1._2.toString())})
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
  }
  def input(name: String) = MyMatrix.input(scala.io.Source.fromFile(name).mkString.split('\n'))
  def testSimple(n: Int): Array[Array[Double]] = {
    import math._
    val A = Array.ofDim[Double](n,n)
    for (i <- 0 until n)
      for (j <- 0 until n) 
        A(i)(j) = Math.sqrt(random*1000)
    A
  }
  def test(n: Int): Array[Array[Double]] = {
    import math._
    val c = (random+2.0)*100.0
    val A = Array.ofDim[Double](n,n)
    for (i <- 0 until n)
      for (j <- 0 until n) 
        A(i)(j) =  c + Math.log1p(i*j) 
    A
  }
}
