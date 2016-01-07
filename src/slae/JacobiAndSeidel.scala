package slae

import structLA._
  
/**
 * @author Ionkin Mikhail
 * This object implements the methods of Jacobi and Seidel
 * [1: 3.2,3.3]

 * [1] -- "С. Ю. Гоголева, Матрицы и вычисления, 2004
 * [2] -- "В. Б. Андреев, Численные методы"
 * [3] -- "J. E. Roman, The QR Algorithm. http://people.inf.ethz.ch/arbenz/ewp/Lnotes/chapter3.pdf"
 * [4] -- "G. Fasshauer  http://www.math.iit.edu/~fass/477577_Chapter_11.pdf"
 */

class JacobiAndSeidel(A: Array[Array[Double]], b: Array[Double]) {      
  val n = A.length
  
  /** It will be stored decision SLAE */ 
  lazy val xsSeidel = {
    val res = new Array[Double](n)
    Array.copy(b, 0, res, 0, n); 
    solveJacobiAndSeidel(res,res)
    res
  }
  
  lazy val xsJacobi = {
    val res = new Array[Double](n)
    Array.copy(b, 0, res, 0, n)
    val xs2 = new Array[Double](n)
    solveJacobiAndSeidel(res,xs2)
    res
  } 
  
  def xsJoS = Function.METHOD_ID match{
    case 1 => xsJacobi
    case 2 => xsSeidel
  }
  
  def methodName = Function.METHOD_ID match{
    case 1 => "Якоби"
    case 2 => "Гаусса-Зейделя"
  }
  
  Counter.countIA+=9; Counter.countMD+=1
  
  /** error of solution */
  var ERROR = 1e-5
  
  lazy val diag = {
    val res = new Array[Double](n)
    for (i <- 0 until n) res(i) = 1.0/A(i)(i)
    res
  }
  
  //  Jacobi
  lazy val BJ  =  {
    val res = new Array[Array[Double]](n);
    for (i <- 0 until n) res(i) = MyVector.multi(A(i), -diag(i))
    for (i <- 0 until n) res(i)(i) +=1
    res
  } 
  // Seidel 
  lazy val BGS = {
    val mat = Array.ofDim[Double](n,n);  
    for (i <- 0 until n; j <- 0 to i) mat(i)(j) = -A(i)(j);
                 
    val R = Array.ofDim[Double](n,n); 
    for (i <- 0 until n; j <- i+1 until n) R(i)(j) = A(i)(j);
                
    val res = MyMatrix.multi(MyMatrix.inverseTriangleDown(mat),R)
    res
  }
  def B = Function.METHOD_ID match{
    case 1 => BJ
    case 2 => BGS
  }
              
  /** q is norm(B) (see [1: p. 21, theor. 2]) */     
  lazy val qJ = Function.normM(BJ,null)
  lazy val qGS = Function.normM(BGS,null)
  def q = Function.METHOD_ID match{
    case 1 => qJ
    case 2 => qGS
  }
  
  /*
   * Необходимое количество итераций можно вычислить заранее!!! 
   * [1: 3.1, th. 2] -- {q = normM(B); dX = normV(x1-x0); e = ERROR;
   * eq: {(q^k)/(1-q)<=e*dX <=> q^k <= e*dX*(1-q), }} 
   * see [1: p. 21, theor. 2]
   */
  lazy val iter  = {
    if (q >= 1.0) throw new java.lang.IllegalArgumentException(
          "Для решения методом "+ methodName +" необходимо слишком много итераций")
    
    val thisdiag = diag.clone()
    for (i <- 0 until n) thisdiag(i)*=b(i)
    val xs2 = MyVector.plus(MyMatrix.multi(B, b), thisdiag)
    val dX = Function.normV(MyVector.minus(xs2, b))
    val res = {
      val m1dX = 1.0/dX;
      var k: Long = 0; var mult = 1.0; 
      if (q<1.0) 
        while (mult>ERROR*(1.0-q)*m1dX) {mult*=q; k+=1}
      else 
        if (mult<ERROR*(1.0-q)*m1dX) {k = 1}
      else k = -1;
      k
    }
    res        
  }
            
  lazy val isConverge = q<1.0
  
  lazy val isStability = {
    Counter.countIA += n+1
    Counter.countC  += n
    
    var flag  = true      // is stable
    for (i <- 0 until n) flag &&= (!MyDouble.isZero(A(i)(i)))
    flag
  }
  
  /** @return solution of SLAE Ax=b (with ERROR) */
  def solveJacobiAndSeidel(xs1: Array[Double], xs2: Array[Double]): Unit = {      
    // [1: 3.2]
    if (!isStability) throw new java.lang.IllegalArgumentException(
      "На диагонали матрицы А есть нуль, решение невозможно")
    if (!isConverge) throw new java.lang.IllegalArgumentException(
      "Решение невозможно -- норма матрицы B не меньше 1.0")
    if (iter == -1) throw new java.lang.IllegalArgumentException(
      "Для решения необходимо слишком много итераций")
        
    Counter.countMD += 2*iter
    Counter.countIA += (n+3)*iter
    Counter.countC  += (n+1)*iter
    Counter.countPM += n*iter
    
    for (j <- 0 until iter.toInt; i <- 0 until n)
      xs2(i) = -(MyVector.multi(A(i),xs1) - A(i)(i)*xs1(i) - b(i))*diag(i)
  }
}