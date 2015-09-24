package slae
/**
 * @author Ionkin Mikhail
 * QR decomposition (also called a QR factorization) 

 * [1] -- "С. Ю. Гоголева, Матрицы и вычисления, 2004
 * [2] -- "В. Б. Андреев, Численные методы"
 * [3] -- "J. E. Roman, The QR Algorithm. http://people.inf.ethz.ch/arbenz/ewp/Lnotes/chapter3.pdf"
 * [4] -- "G. Fasshauer  http://www.math.iit.edu/~fass/477577_Chapter_11.pdf"
 */
class QR(A: Array[Array[Double]], b: Array[Double]) {
      import LR1._
      import LR1.MyDouble
      import math._
      val n = A.length
      Counter.countIA+=2;
      var ERROR = 1e-3
      var ITER = 10
      /** 
       *  Householder reflection 
       *  [1: 1.3--Метод отражений]
       *  */
      def getQR_HR(A: Array[Array[Double]]): (Array[Array[Double]],Array[Array[Double]]) = {
        val n = A.length
        Counter.countIA+=(n*n)<<2; Counter.countC+=(n*n)<<2;
        
        var curH1: Array[Array[Double]] = null; val curH2 = Array.ofDim[Double](n,n)
        var curR = Array.ofDim[Double](n,n); MyMatrix.copy(A,curR)
        
        solve()
        
        def solve(){
          val normID = ISLAE.NORM_ID
          ISLAE.NORM_ID = 2
          Counter.countC+=n-1;
          for (k <- 0 until n-1){            
            Counter.countIA+=((n-k)<<1); Counter.countC+=((n-k)<<1); Counter.countMD+=1; Counter.countIA+=4;
            
            val v = new Array[Double](n-k); for (i <- k until n) v(i-k) = curR(i)(k); v(0) += MyDouble.sign(v(0))*ISLAE.normV(v)     // p = v; [1, 1.3 -- eq 1.14]
            curH1 = getHouseholdersMatrix(v)
            setH2()

            curR = MyMatrix.multi(curH2,curR)            
          }
          ISLAE.NORM_ID = normID
          /**
           * Householders Matrix
           */
          def getHouseholdersMatrix(p: Array[Double]): Array[Array[Double]] = {
            Counter.countMD+=1
            val matrix = MyMatrix.multi(MyVector.multiGM(p,p),-2.0/MyVector.multi(p,p))
            Counter.countIA += matrix.length; Counter.countC += matrix.length; Counter.countPM+=matrix.length;
            for (i <- 0 until matrix.length) matrix(i)(i) += 1
            matrix
          }
          def setH2(){
            val i2 = n - curH1.length
            Counter.countC += n*n; Counter.countIA+= n*n;
            for (i <- 0 until i2) for (j <- 0 until i2) curH2(i)(j) = if (i == j) 1.0 else 0.0
            for (i <- 0 until i2) for (j <- i2 until n) curH2(i)(j) = 0.0
            for (j <- 0 until i2) for (i <- i2 until n) curH2(i)(j) = 0.0
            for (i <- i2 until n) for (j <- i2 until n) curH2(i)(j) = curH1(i-i2)(j-i2)
          }
        }
        val Q = MyMatrix.multi(A, MyMatrix.inverseTriangleUp(curR));//LU.getAndSetInverseMatrix(curR))
        (Q,curR)
      }
      // Givens rotations
      def getQR_GR(A: Array[Array[Double]] ): (Array[Array[Double]],Array[Array[Double]]) = { 
        val n = A.length
        Counter.countIA+=n*n*3; Counter.countC+=n*n*3
        var curR = Array.ofDim[Double](n,n); MyMatrix.copy(A, curR)        
        for (i <- 0 until n){
          Counter.countC+=(i+1)<<1
          for (j <- 0 until i) 
            if (!MyDouble.isZero(curR(i)(j))) {
              Counter.countIA+=(n<<1)+3;
              val colm = new Array[Double](n); for (k <- 0 until n) colm(k)=curR(k)(j)
              val (c,s) = getGivensMatrixCS(colm,i,j)
              curR = multiGR(curR,c,s,i,j)
            }
        }
        /** 
         *  (c,s)
         **/
        def getGivensMatrixCS(v: Array[Double], i: Int, j: Int): (Double, Double) = {
          val n = v.length
          Counter.countMD+=5
          val m1tay = 1.0/sqrt(v(i)*v(i)+v(j)*v(j))
          (v(j)*m1tay, v(i)*m1tay)
        }
        /**
        * 1 0  0 0    11 12 13 14   11      12      13      14 
        * 0 c  s 0  * 21 22 23 24 = c21+s31 c22+s32 c23+s33 c24+s34 
        * 0 -s c 0    31 32 33 34  -s21+c31-s22+c32
        * 0 0  0 1    41 42 43 44   41      42      43      44
        **/
        def multiGR(B: Array[Array[Double]], c: Double, s: Double,  i00: Int, j00: Int): Array[Array[Double]] = {
          val n = A.length
          val res = Array.ofDim[Double](n,n)
          for (i <- 0 until n) if (i!=i00 && i!=j00) for (j <- 0 until n) res(i)(j) = B(i)(j)
          var i0 = min(i00,j00); var j0 = max(i00,j00); 
          for (j <- 0 until n) res(i0)(j) = c*B(i0)(j)+s*B(j0)(j)
          for (j <- 0 until n) res(j0)(j) = -s*B(i0)(j)+c*B(j0)(j)
          res
        }
        val Q = MyMatrix.multi(A, MyMatrix.inverseTriangleUp(curR))
        (Q,curR)
      }
      val Ak = getAk() 
      /**
       * Basic QR-algorithm
       * [3, p. 52]
       */
      private def getAk(): Array[Array[Double]] = {  //_WithShifts(): Array[Array[Double]] = {  
        Counter.countIA+=n*n+8
        val updA2 = Array.ofDim[Double](n,n); MyMatrix.copy(A, updA2); MyMatrix.transp(updA2); 
        val updA = MyMatrix.multi(A,updA2); //   Array.ofDim[Double](n,n); MyMatrix.copy(A, updA); //
        // var Ann = updA(n-1)(n-1); for (i <- 0 until n) updA(i)(i)-=Ann
        var QR = getQR_HR(updA)
        var curA = MyMatrix.multi(QR._2,QR._1); // for (i <- 0 until n) curA(i)(i)+=Ann
        var norm2 = ISLAE.normM(curA,null)
        Counter.countMD+=1
        var norm1 = norm2+ERROR*10+1.0
        var iter = ITER
        while (abs(norm2-norm1)>ERROR && iter>0){
          Counter.countIA+=6; Counter.countC+=1
          // var Ann = updA(n-1)(n-1); for (i <- 0 until n) updA(i)(i)-=Ann
          QR = getQR_GR(curA)
          curA = MyMatrix.multi(QR._2,QR._1); // for (i <- 0 until n) curA(i)(i)+=Ann
          norm1 = norm2
          norm2 = ISLAE.normM(curA,null)
          iter-=1
        }
        if (iter == 0) throw new java.lang.IllegalArgumentException(
            "Последовательность {Ak}, k=0:infinity не сходится по норме на "+ITER.toInt+" итерациях (с точностью до ERROR = "+ERROR+").\n"+
            "Предполагаемые собственные числа: "+
                MyVector.toString(
                  {
                    val diag = new Array[Double](n)
                    for (i <- 0 until n) diag(i) = curA(i)(i)
                    diag
                  }))
        curA
      }
      /** 
       *  bad work!!!
       *  QR algorithm with shifts
       *  [3, p. 63: Algorithm 3.4 The Hessenberg QR algorithm with Rayleigh quotient shift 
       *  //[3, p. 70 Algorithm 3.5 The Francis double step QR algorithm]
       *  //H_ij ??? 
       *  //P ???
       *  
       *   */
      def AkWithShifts(): Array[Array[Double]] = { 
        val updA2 = Array.ofDim[Double](n,n); MyMatrix.copy(A, updA2); MyMatrix.transp(updA2); 
        val updA = MyMatrix.multi(A,updA2);
        var curA = Array.ofDim[Double](n,n); MyMatrix.copy(updA, curA);
        /** [3, p. 63: Algorithm 3.4 The Hessenberg QR algorithm with Rayleigh quotient shift */
        { 
          var QR: (Array[Array[Double]],Array[Array[Double]]) = null
          var sigmaK = 0.0
          var m = n-1; 
          QR = getQR_HR(curA); 
          curA = MyMatrix.multi(QR._2,QR._1); 
          while (m>=1){ 
            do {
              sigmaK = curA(m)(m);                 for (i <- 0 until n) curA(i)(i) -= sigmaK
              QR = getQR_GR(curA); 
              curA = MyMatrix.multi(QR._2,QR._1);  for (i <- 0 until n) curA(i)(i) += sigmaK
            }
            while (abs(curA(m)(m-1))>1e-1);
            m -= 1
          }
        }
        curA
      }
      def getSolve(QR: (Array[Array[Double]], Array[Array[Double]])): Array[Double] = {
        Counter.countIA+=4+(n<<2); Counter.countPM+=(n<<1)
        val Q = QR._1
        val R = QR._2
        val newB = new Array[Double](n)
        for (k <- 0 until n) newB(k) = MyVector.multi(b,Q,k)
        val xs = new Array[Double](n)
        // index
        var i = n-1
        while (i>=0) {
          Counter.countMD+=2
          xs(i) = (1.0/R(i)(i))*(newB(i) - MyVector.multi(R(i),xs,i+1,n))
          i-=1;
        }
        xs
      }
      val eigenvalues = {          
        Counter.countC+=n; Counter.countIA+=n;
        val diag = new Array[Double](n)
        for (i <- 0 until n) diag(i) = Ak(i)(i)
        diag
      }
}
