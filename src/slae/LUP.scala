package slae

/**
 * @author Ionkin Mikhail
 * LU factorization with Partial Pivoting 
 * [1 : 1.2], [2: 4.4] 

 * [1] -- "С. Ю. Гоголева, Матрицы и вычисления, 2004
 * [2] -- "В. Б. Андреев, Численные методы"
 * [3] -- "J. E. Roman, The QR Algorithm. http://people.inf.ethz.ch/arbenz/ewp/Lnotes/chapter3.pdf"
 * [4] -- "G. Fasshauer  http://www.math.iit.edu/~fass/477577_Chapter_11.pdf"
 */
class LUP (A: Array[Array[Double]], b: Array[Double]) {
  val n = A.length  
  lazy val LU = Function.LU.Decomposition.LU(A, true, q);
  lazy val L = LU._1
  lazy val U = LU._2
  // the resulting vector permutations
  lazy val q = {
    val res = new Array[Int](n); for (i <- 0 until n) res(i) = i; res
  }
  
  /** solution of SLAE Ax = b */
  lazy val xs = Function.LU.Solution.solve(L, U, b, q)
  
  /**
  * A = LU; 
  * [1, eqs. 1.8-1.10]
  * @return A^(-1)
  */ 
  lazy val inverse = Function.LU.Info.inverse(L, U)
  lazy val cond    = Function.LU.Info.cond(A, L, U)
  lazy val det     = Function.LU.Info.det(U)
}