package slae
<<<<<<< HEAD

=======
>>>>>>> origin/master
/**
 * @author Ionkin Mikhail
 * LU factorization with Partial Pivoting 
 * [1 : 1.2], [2: 4.4] 

 * [1] -- "С. Ю. Гоголева, Матрицы и вычисления, 2004
 * [2] -- "В. Б. Андреев, Численные методы"
 * [3] -- "J. E. Roman, The QR Algorithm. http://people.inf.ethz.ch/arbenz/ewp/Lnotes/chapter3.pdf"
 * [4] -- "G. Fasshauer  http://www.math.iit.edu/~fass/477577_Chapter_11.pdf"
 */
<<<<<<< HEAD
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
=======
class LUP(A: Array[Array[Double]], b: Array[Double]) {
      import LR1._
      val n = A.length
      // copy of A
      private def matr = {
        Counter.countIA+=(n*n)<<1; Counter.countC+=(n*n)<<1; 
        val matr = Array.ofDim[Double](n,n); MyMatrix.copy(A, matr); matr
      }
      Counter.countC+=n; Counter.countIA+=n+7;
      // the resulting vector permutations
      private val q = new Array[Int](n); for (i <- 0 until n) q(i) = i
      /** A =: LU*/
      val LU: (Array[Array[Double]],Array[Array[Double]]) = ISLAE.getLU(matr, true, q)
      /** det(A) */
      def det = MyMatrix.detDiagMatr(LU._2)
      /** solution of SLAE Ax = b */
      def xs = {
        Counter.countIA+=1        
        // необходимо найти перестановку, обратную к q, и преобразовать решение:
        val x = ISLAE.solve(LU._1,LU._2,b);
        MyVector.swapInverse(x,q); x
      }
}
>>>>>>> origin/master
