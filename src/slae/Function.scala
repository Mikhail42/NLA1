package slae

import structLA._
import math._

/**
 * @author Ionkin Mikhail
 * some static methods for SLAE

 * [1] -- "С. Ю. Гоголева, Матрицы и вычисления, 2004
 * [2] -- "В. Б. Андреев, Численные методы"
 * [3] -- "J. E. Roman, The QR Algorithm. http://people.inf.ethz.ch/arbenz/ewp/Lnotes/chapter3.pdf"
 * [4] -- "G. Fasshauer  http://www.math.iit.edu/~fass/477577_Chapter_11.pdf"
 */
object Function {
    
     /** Jacobi or Seidel 
     *  1 - Jacobi method
     *  2 - Seidel method
     *  */
    var METHOD_ID = 2

    /** 
     *  p of p-norms: p = NORM_ID
     *  if (p==3) => (p = Infinity)
     *  */
    var NORM_ID = 3
    
    var QR_ERROR = 1e-3
    
    var QR_ITER = 10
    
    /**
     * vector norm
     */
    def normV(a: Array[Double]) = NORM_ID match{
      case 1 => MyVector.norm1(a)
      case 2 => MyVector.norm2(a)
      case 3 => MyVector.normCubic(a)
    }
    
    /**
     * matrix norm
     */
    def normM(matr: Array[Array[Double]], eigenvalues: Array[Double]) = NORM_ID match{
      case 1 => MyMatrix.normMC(matr)
      case 2 => {
        if (eigenvalues == null) throw new java.lang.IllegalArgumentException(
            """Неверный аргумент. На данном этапе невозможно (эффективно) посчитать собственные числа (проверьте это). 
              Пожалуйста, выберите другую норму или по-другому реализуйте методы""");
        else MyMatrix.normSpectr(eigenvalues)
      }
      case 3 => MyMatrix.normMS(matr)
    }
    
  /**
  * norm(a.-b)/norm(a)
  */
  def error(a: Array[Double], b: Array[Double]): Double = {
    Counter.countMD+=1; Counter.countIA+=1
    val c = MyVector.minus(a, b)
    normV(c)/normV(a)
  }   
    
  object LLT{
    object Decomposition {
      def getLLT(matr: Array[Array[Double]]): Array[Array[Double]]= {
        val n = matr.length
        import math._
        val L = Array.ofDim[Double](n,n)
        for (i <- 0 until n)
          for (j <- 0 to i)
            if (i==j) L(i)(i) = sqrt(matr(i)(i)-MyVector.multi(L(i),L(i),0,i))
            else L(i)(j) = (1.0/L(j)(j))*(matr(i)(j)-MyVector.multi(L(i),L(j),0,j))
        L
      }
    }
  }
  
  object LU{
    
    object Info {
      
      def det(U: Array[Array[Double]]): Double = {
          if (U == null) throw new java.lang.NullPointerException(
              "Для выполнения данной операции необходимо выполнить LU-разложение")
          MyMatrix.detDiagMatr(U)
      }
      
      def cond(A:Array[Array[Double]], L: Array[Array[Double]], U: Array[Array[Double]]): Double = {
          if (L == null || U == null) throw new java.lang.NullPointerException(
              "Для выполнения данной операции необходимо выполнить LU-разложение")
          (MyMatrix.normMS(A)*MyMatrix.normMS(inverse(L,U)))
      }
    
      /**
       * A = LU; 
       * [1, eqs. 1.8-1.10]
       * @return A^(-1)
      */
      def inverse(L:Array[Array[Double]], U: Array[Array[Double]]): Array[Array[Double]] = {
        if (L == null || U == null) throw new java.lang.NullPointerException("Для выполнения данной операции необходимо выполнить LU-разложение")
        val n = L.length
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
    }
    
    object Solution {
      
      /**
      * [1], p. 7, eq. (1.5)
      */
      def solveLUP(L: Array[Array[Double]], U: Array[Array[Double]], b: Array[Double]): Array[Double] = {
        if (L == null || U == null) throw new java.lang.NullPointerException("Для выполнения данной операции необходимо выполнить LU-разложение")
        val n = L.length
        Counter.countMD+=n; Counter.countIA+=(n<<2); Counter.countPM+=(n+n); Counter.countC+=(n+n)
        val ys = new Array[Double](n)
        import structLA.MyVector.multi
        for (i <- 0 until n) ys(i) = b(i) - multi(L(i),ys,0,i)
        val xs = new Array[Double](n)
        var i = n-1
        while (i>=0) {
          xs(i) = (1.0/U(i)(i))*(ys(i) - multi(U(i),xs,i+1,n))
          i-=1
        }
        xs
      }
      
      /** solution of SLAE Ax = b */
      def solve(L:Array[Array[Double]], U: Array[Array[Double]], b:Array[Double] , q: Array[Int]) =  {
        if (L == null || U == null) throw new java.lang.NullPointerException(
            "Для выполнения данной операции необходимо выполнить LU-разложение")
        Counter.countIA+=1        
        // необходимо найти перестановку, обратную к q, и преобразовать решение:
        val x = solveLUP(L, U, b);
        MyVector.swapInverse(x,q); 
        x
      }  
    }
    
    object Decomposition {
      /**
       * LU-decomposition. 
       * [1 : 1.1,1.2], [2: 1.2, 4.4] 
       * @param dwpFlag -- (true => decomposition with pivoting, false => simple decomposition)
       * @q -- the resulting vector permutations (if dwpFlag == true). 
       * @return (L,U). q is UPDATE
       */
      @throws(classOf[java.lang.IllegalArgumentException])
      def LU(matr: Array[Array[Double]], dwpFlag: Boolean, q: Array[Int]): (Array[Array[Double]],Array[Array[Double]]) = {       
        val n = matr.length
        val U = Array.ofDim[Double](n,n)
        val L = Array.ofDim[Double](n,n)
        Counter.countIA+=n*n; Counter.countC+=n*n+n; Counter.countPM+=n*n
        //ДлConsole.println(MyMatrix.toString(matr))
        /**
         *  [2]: 
         *   стр. 11: Можно чередовать вычисление строк U(k)(:) и столбцов L(:)(k)
         *   стр. 42: Устойчивость не гарантирована 
         */
        for (k <- 0 until n){
          for (j <- k until n) U(k)(j) = matr(k)(j) - MyVector.multi(L(k),U,j,0,k);
          if (dwpFlag){
            // индекс элемента с наибольшим модулем
            val ind = MyVector.indexOfMaxAbs(matr(k),k) // val ind = MyVector.indexOfMaxAbs(U(k),k)
           // Console.println(k+" "+ind)
            MyMatrix.swapColumn(matr,k,ind);  MyMatrix.swapColumn(U,k,ind)
            Counter.countC+=1
            if (k!=ind) {Counter.countIA+=3;  val g = q(k);  q(k)=q(ind);  q(ind) = g;}
          }
          if (MyDouble.isZero(U(k)(k))) throw new java.lang.IllegalArgumentException(
              "U("+k+")("+k+")="+U(k)(k)+". \n\t\t\t LU-разложение невозможно ввиду необходимости деления на сверхмалое число")
          Counter.countMD+=(n-k); Counter.countIA+=(n-k+1); Counter.countPM+=(n-k); Counter.countC+=(n-k)
          val m1Ukk = 1.0/U(k)(k); L(k)(k) = 1.0
          for (i <- k+1 until n) L(i)(k) = (matr(i)(k) - MyVector.multi(L(i),U,k,0,k))*m1Ukk
        }      
        (L,U)
      }
    }

  }
  
  object QR {  
    
    var ERROR = Function.QR_ERROR
    var ITER = Function.QR_ITER
    
    import Eignvalues._
    import Decomposition._
    
    object Solution {
      
      def solve(QR: (Array[Array[Double]], Array[Array[Double]]), b: Array[Double]): Array[Double] = {
        val n = QR._1.length
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
    }
    
    object Eignvalues {
      
      /**
       * Basic QR-algorithm
       * [3, p. 52]
       */
      private def getAk(A: Array[Array[Double]]): Array[Array[Double]] = {  //_WithShifts(): Array[Array[Double]] = {  
        val n = A.length
        Counter.countIA+=n*n+8
        val updA2 = Array.ofDim[Double](n,n); MyMatrix.copy(A, updA2); MyMatrix.transp(updA2); 
        val updA = MyMatrix.multi(A,updA2); //   Array.ofDim[Double](n,n); MyMatrix.copy(A, updA); //
        // var Ann = updA(n-1)(n-1); for (i <- 0 until n) updA(i)(i)-=Ann
        var QR = getQR_HR(updA)
        var curA = MyMatrix.multi(QR._2,QR._1); // for (i <- 0 until n) curA(i)(i)+=Ann
        var norm2 = Function.normM(curA,null)
        Counter.countMD+=1
        var norm1 = norm2+ERROR*10+1.0
        var iter = ITER
        while (abs(norm2-norm1)>ERROR && iter>0){
          Counter.countIA+=6; Counter.countC+=1
          // var Ann = updA(n-1)(n-1); for (i <- 0 until n) updA(i)(i)-=Ann
          QR = getQR_GR(curA)
          curA = MyMatrix.multi(QR._2,QR._1); // for (i <- 0 until n) curA(i)(i)+=Ann
          norm1 = norm2
          norm2 = Function.normM(curA,null)
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
      private def AkWithShifts(A: Array[Array[Double]]): Array[Array[Double]] = { 
        val n = A.length
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
      
      def eigenvalues(A: Array[Array[Double]]) = {  
        val n = A.length
        val Ak = getAk(A) 
        Counter.countC +=n ; Counter.countIA += n+1;
        val diag = new Array[Double](n)
        for (i <- 0 until n) diag(i) = Ak(i)(i)
        diag
      }
    }
    
    object Decomposition {
      /** 
      *  Householder reflection 
      *  [1: 1.3--Метод отражений]
      *  */
      def getQR_HR(A: Array[Array[Double]]): (Array[Array[Double]],Array[Array[Double]]) = {
        val n = A.length
        Counter.countIA+=(n*n)<<2; Counter.countC+=(n*n)<<2;
           
        var curH1: Array[Array[Double]] = null; val curH2 = Array.ofDim[Double](n,n)
        var curR = Array.ofDim[Double](n,n); MyMatrix.copy(A,curR)
            
        val normID = Function.NORM_ID
        Function.NORM_ID = 2
        Counter.countC+=n-1;
        for (k <- 0 until n-1){            
          Counter.countIA+=((n-k)<<1); Counter.countC+=((n-k)<<1); Counter.countMD+=1; Counter.countIA+=4;
            
          val v = new Array[Double](n-k); 
          for (i <- k until n) v(i-k) = curR(i)(k); 
          v(0) += MyDouble.sign(v(0))*Function.normV(v)     // p = v; [1, 1.3 -- eq 1.14]
                
          curH1 = getHouseholdersMatrix(v)
          setH2()
    
          curR = MyMatrix.multi(curH2,curR)            
        }
        Function.NORM_ID = normID
        
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
          for (i <- 0 until i2; j <- 0 until i2) curH2(i)(j) = if (i == j) 1.0 else 0.0 
          for (i <- 0 until i2) for (j <- i2 until n) curH2(i)(j) = 0.0
          for (j <- 0 until i2; i <- i2 until n) curH2(i)(j) = 0.0
          for (i <- i2 until n; j <- i2 until n) curH2(i)(j) = curH1(i-i2)(j-i2)
        }
        
        val Q = MyMatrix.multi(A, MyMatrix.inverseTriangleUp(curR));
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
          for (j <- 0 until n) res(i0)(j) =  c*B(i0)(j)+s*B(j0)(j)
          for (j <- 0 until n) res(j0)(j) = -s*B(i0)(j)+c*B(j0)(j)
          res
        }
        val Q = MyMatrix.multi(A, MyMatrix.inverseTriangleUp(curR))
        (Q,curR)
      }
    }
  }
}