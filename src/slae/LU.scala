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
import LR1._
class LU(A: Array[Array[Double]], b: Array[Double]) {
      val n = A.length
      Counter.countIA+=5
      @throws(classOf[java.lang.IllegalArgumentException])
      private val matrixsLU = ISLAE.getLU(A, false, null)
      val L = matrixsLU._1
      val U = matrixsLU._2
      val xs = ISLAE.solve(L,U,b)
      /**
       * A = LU; 
       * [1, eqs. 1.8-1.10]
       * @param U
       * @return A^(-1)
      */
      def getAndSetInverseMatrix(L:Array[Array[Double]], U: Array[Array[Double]]): Array[Array[Double]] = {
        import MyVector.multi
        Counter.countIA+=n*n+1
        val mat = Array.ofDim[Double](n,n);
        /** [1, eqs. 1.8-1.10] */
        def setMat(i: Int, j: Int) : Unit = {
          Counter.countIA+=1; Counter.countMD+=2 // max of probability
          mat(i)(j) = 
            if (i == j) (1.0 - multi(U(i),mat,j,j+1,n))/U(j)(j)
            else if (i < j) -multi(U(i),mat,j,i+1,n)/U(i)(i)
            else -multi(U(i),L,j,j+1,n)
        }
        Counter.countIA+=n*n; Counter.countC==n*n; 
        var i = n-1; while (i>=0){var j = n-1; while (j>=0) {setMat(i,j); j-=1;}; i-=1;}

        mat
      }
      Counter.countMD+=1; Counter.countIA+=2;
      def getAndSetCond(): Double = (MyMatrix.normMS(A)*MyMatrix.normMS(getAndSetInverseMatrix(L,U)))
      def getAndSetDet(): Double = MyMatrix.detDiagMatr(U)
}
