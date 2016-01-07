package slae
/**
 * @author Ionkin Mikhail
 * uses LU-decomposition
 * [1: 1.1], [2: 1.2] 

 * [1] -- "С. Ю. Гоголева, Матрицы и вычисления, 2004
 * [2] -- "В. Б. Андреев, Численные методы"
 * [3] -- "J. E. Roman, The QR Algorithm. http://people.inf.ethz.ch/arbenz/ewp/Lnotes/chapter3.pdf"
 * [4] -- "G. Fasshauer  http://www.math.iit.edu/~fass/477577_Chapter_11.pdf"
 */
class LU(A: Array[Array[Double]], b:  Array[Double]){ 
  val n = A.length  
  lazy val matrixsLU = Function.LU.Decomposition.LU(A, false, null)
  lazy val L = matrixsLU._1
  lazy val U = matrixsLU._2
  
  /** solution of SLAE Ax = b */
  lazy val xs = Function.LU.Solution.solveLUP(L, U, b)
  
  /**
  * A = LU; 
  * [1, eqs. 1.8-1.10]
  * @return A^(-1)
  */ 
  lazy val inverse = Function.LU.Info.inverse(L, U)
  lazy val cond    = Function.LU.Info.cond(A, L, U)
  lazy val det     = Function.LU.Info.det(U)
}