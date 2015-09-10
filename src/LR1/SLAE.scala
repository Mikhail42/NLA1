package LR1;
/**
 * @author Ionkin Mikhail
 */
class SLAE(A : Array[Array[Double]], x: Array[Double]){
    import math._
    import MyDouble._
    // for p-norms, where p = NORM_ID
    var NORM_ID = 3
    // Jacobi or Seidel
    var METHOD_ID = 1
    val n = x.length
    val b = MyMatrix.multi(A,x)
    def normV(a: Array[Double]) = NORM_ID match{
      case 1 => MyVector.norm1(a)
      case 2 => MyVector.norm2(a)
      case 3 => MyVector.normCubic(a)
    }
    def normM(matr: Array[Array[Double]], eigenvalues: Array[Double]) = NORM_ID match{
      case 1 => MyMatrix.normMC(matr)
      case 2 => {
        if (eigenvalues == null) throw new java.lang.IllegalArgumentException(
            """Неверный аргумент. На данном этапе невозможно (эффективно) посчитать собственные числа (проверьте это). 
              Пожалуйста, выберите другую норму""");
        else MyMatrix.normSpectr(eigenvalues)
      }
      case 3 => MyMatrix.normMS(matr)
    }
    /**
     * the relative error SLAE solutions using standards normCubic (see object MyVector)
    */
    def error(a: Array[Double], b: Array[Double]): Double = {
      val c = MyVector.minus(a, b)
      normV(c)/normV(a)
    }
    @throws(classOf[java.lang.IllegalArgumentException])
    private def getLU(matr: Array[Array[Double]], dwpFlag: Boolean, q: Array[Int]): (Array[Array[Double]],Array[Array[Double]]) = {       
      import MyMatrix._
      val U = Array.ofDim[Double](n,n)
      val L = Array.ofDim[Double](n,n)
      /**
       *  В. Б. Андреев, Численные методы
       *   стр. 11: Можно чередовать вычисление строк U(k)(:) и столбцов L(:)(k)
       *   стр. 42: Устойчивость не гарантирована 
       */
      for (k <- 0 until n){
        for (j <- 0 until n) U(k)(j) = matr(k)(j) - MyVector.multi(L(k),U,j,0,k);
        if (dwpFlag){
          // индекс элемента с наибольшим модулем
          val ind = MyVector.indexOfMaxAbs(U(k),k)
          swapColumn(matr,k,ind); swapColumn(U,k,ind); 
          if (k!=ind) {val g = q(k); q(k)=q(ind); q(ind) = g;}
        }
        if (isZero(U(k)(k))) throw new java.lang.IllegalArgumentException("LU-разложение невозможно ввиду необходимости деления на сверхмалое число")
        Counter.countOpsMD+=(n-k)
        val m1Ukk = 1.0/U(k)(k); L(k)(k) = 1.0;
        for (i <- k+1 until n) L(i)(k) = (matr(i)(k) - MyVector.multi(L(i),U,k,0,k))*m1Ukk; 
      }
      (L,U)
    }
    private def solve(L: Array[Array[Double]], U: Array[Array[Double]]): Array[Double] = {
      val ys = new Array[Double](n)
      import MyVector.multi
      for (i <- 0 until n) ys(i) = b(i) - multi(L(i),ys,0,i)
      val xs = new Array[Double](n)
      var i = n-1
      while (i>=0) {
        xs(i) = (1.0/U(i)(i))*(ys(i) - multi(U(i),xs,i+1,n))
        i-=1
      }
      xs
    }
    /**
     * uses LU-decomposition
     */
    object LU{  
      @throws(classOf[java.lang.IllegalArgumentException])
      private val matrixsLU = getLU(A, false, null)
      val L = matrixsLU._1
      val U = matrixsLU._2
      val xs = solve(L,U)

      val inverseMatrix = {
        import MyVector.multi
        val mat = Array.ofDim[Double](n,n);
        def setMat(i: Int, j: Int) : Unit = {
          Counter.countOpsMD+=2 // max of probability
          mat(i)(j) = 
            if (i == j) (1.0/U(j)(j))*(1.0 - multi(U(i),mat,j,j+1,n))
            else if (i < j) -(1.0/U(i)(i))*multi(U(i),mat,j,i+1,n)
            else -multi(U(i),mat,j,j+1,n)
        }
        var i = n-1; while (i>=0){var j = n-1; while (j>=0) {setMat(i,j); j-=1;}; i-=1;}
        mat
      }
      Counter.countOpsMD+=1
      val cond = MyMatrix.normMS(A)*MyMatrix.normMS(inverseMatrix)
      val det = MyMatrix.detDiagMatr(U)
    }
    /**
     * LU - decomposition with pivoting
     */
    object LU_DWP{
      private val matr = {
        val matr = Array.ofDim[Double](n,n);    
        for (i <- 0 until n)  Array.copy(A(i), 0, matr(i), 0, n)
        matr
      }
      val q = new Array[Int](n); for (i <- 0 until n) q(i) = i
      private val matrixsLU =  getLU(matr, true, q)
      val L = matrixsLU._1
      val U = matrixsLU._2
      val det = MyMatrix.detDiagMatr(U)
      val xs = {
        val x = solve(L,U)
        // необходимо найти перестановку, обратную к q, и преобразовать решение:
        MyVector.swapInverse(x,q)
        x        
      }
    }
    object JacobiAndSeidel{
      var ERROR = 1e-5
      var ITER = 1e4.toInt
      private val xs = new Array[Double](n); Array.copy(b, 0, xs, 0, n)
      private var norm0 = normV(xs)
      Counter.countOpsMD+=1
      private var norm = norm0+10.0*ERROR+0.01;
      private val m1diagA = new Array[Double](n); {Counter.countOpsMD+=n; for (i <- 0 until n) m1diagA(i) = 1.0/A(i)(i);}
      private val isStability = {var flag  = true; for (i <- 0 until n) flag&&=(MyDouble.compareDouble(A(i)(i),0.0) != 0); flag}
      private def solve(xs1: Array[Double], xs2: Array[Double]): Unit = {
        val res = new Array[Double](n)
        if (!isStability) throw new java.lang.IllegalArgumentException("На диагонали матрицы А есть нуль, решение невозможно")
        var iter = ITER
        while (abs(norm-norm0)>ERROR && iter>0){
          Counter.countOpsMD+=2
          for (i <- 0 until n) xs2(i) = -(MyVector.multi(A(i),xs1) - A(i)(i)*xs1(i) - b(i))*m1diagA(i)
          norm0 = norm
          norm = normV(xs2)
          iter-=1
        }
        if (iter == 0) throw new java.lang.IllegalArgumentException(
            "\n   Решение по методу "+ (if (METHOD_ID == 1) "Якоби" else "Зейдаля")+
            " не сходится по норме на "+ITER.toInt+" итерациях (с точностью до ERROR = "+ERROR+
            ").\nПредполагаемый ответ: "+MyVector.toString(xs)+
            ". Вы можете изменить число итераций или точноcть, \nобратившись к полям ITER и/или ERROR объекта JacobiAndSeidel.")
      }
      //Philipp Ludwig von Seidel
      val xsSeidel = {solve(xs,xs); xs}
      val xsJacobi = {val xs2 = new Array[Double](n); solve(xs,xs2); xs}
    }
    object QR{      
      var ERROR = 1e-3
      var ITER = 1e4.toInt
      // Householder reflection
      // val QR_HR  = getQR_HR(A)
      def getQR_HR(A: Array[Array[Double]]): (Array[Array[Double]],Array[Array[Double]]) = {
        val n = A.length
        var curH1: Array[Array[Double]] = null; 
        val curH2 = Array.ofDim[Double](n,n)
        var curQ = Array.ofDim[Double](n,n)
        var curR = Array.ofDim[Double](n,n)
        MyMatrix.copy(A,curR)
        solve()
        def solve(){
          val normID = NORM_ID
          NORM_ID = 2
          for (k <- 0 until n-1){
            val v = new Array[Double](n-k)
            for (i <- k until n) v(i-k) = curR(i)(k)
            Counter.countOpsMD+=1
            v(0) += sign(v(0))*normV(v)     // p = v
            curH1 = getHouseholdersMatrix(v)
            setH2()
            if (k == 0) MyMatrix.copy(curH2,curQ)
            else curQ = MyMatrix.multi(curH2,curQ)
            curR = MyMatrix.multi(curQ,A)
          }
          NORM_ID = normID;
          /**
           * Householders Matrix
           */
          def getHouseholdersMatrix(p: Array[Double]): Array[Array[Double]] = {
            Counter.countOpsMD+=1
            val matrix = MyMatrix.multi(MyVector.multiGM(p,p),-2.0/MyVector.multi(p,p))
            for (i <- 0 until matrix.length) matrix(i)(i) += 1
            matrix
          }
          def setH2(){
            val i2 = n - curH1.length
            for (i <- 0 until i2) for (j <- 0 until i2) curH2(i)(j) = if (i == j) 1.0 else 0.0
            for (i <- 0 until i2) for (j <- i2 until n) curH2(i)(j) = 0.0
            for (j <- 0 until i2) for (i <- i2 until n) curH2(i)(j) = 0.0
            for (i <- i2 until n) for (j <- i2 until n) curH2(i)(j) = curH1(i-i2)(j-i2)
          }
        }
        MyMatrix.transp(curQ)
        (curQ,curR)
      }
      // Givens rotations
      def getQR_GR(A: Array[Array[Double]] ): (Array[Array[Double]],Array[Array[Double]]) = {  
        def getGivensMatrix(v: Array[Double], i: Int, j: Int): Array[Array[Double]] = {
          val n = v.length
          Counter.countOpsMD+=5
          val m1tay = 1.0/sqrt(v(i)*v(i)+v(j)*v(j))
          val c = v(j)*m1tay
          val s = v(i)*m1tay
          val G = Array.ofDim[Double](n,n)
          for (k <- 0 until n) if (k!=i && k!=j) G(k)(k) = 1.0;
          G(i)(i) = c; G(j)(j) = c; G(max(i,j))(min(i,j)) = -s; G(min(i,j))(max(i,j)) = s;
          G
        }
        import MyMatrix._
        val n = A.length
        var curQ = Array.ofDim[Double](n,n); 
        var curR = Array.ofDim[Double](n,n)
        var flagQ = false
        MyMatrix.copy(A, curR)
        for (i <- 0 until n){
          for (j <- 0 until i) 
            if (compareDouble(curR(i)(j),0.0) != 0) {
              val colm = new Array[Double](n)
              for (k <- 0 until n) colm(k)=curR(k)(j)
              val G = getGivensMatrix(colm,i,j)
              if (!flagQ){
                MyMatrix.copy(G,curQ)
                flagQ = true          
              } else curQ = multi(G,curQ)
              curR = multi(G,curR)
            }
        }
        if (!flagQ) for (i <- 0 until n) curQ(i)(i) = 1.0
        transp(curQ)
        (curQ,curR)
      }
      val Ak = getAk_WithShifts()
      /**
       * algorithm with shifts
       */
      private def getAk_WithShifts(): Array[Array[Double]] = {  
        val updA = A // Array.ofDim[Double](n,n); MyMatrix.copy(A, updA)
        // var Ann = updA(n-1)(n-1); for (i <- 0 until n) updA(i)(i)-=Ann
        var QR = getQR_HR(updA)
        var curA = MyMatrix.multi(QR._2,QR._1); // for (i <- 0 until n) curA(i)(i)+=Ann
        var norm2 = normM(curA,null)
        Counter.countOpsMD+=1
        var norm1 = norm2+ERROR*10+1.0
        var iter = ITER
        while (abs(norm2-norm1)>ERROR && iter>0){
          var Ann = updA(n-1)(n-1); // for (i <- 0 until n) updA(i)(i)-=Ann
          QR = getQR_GR(curA)
          curA = MyMatrix.multi(QR._2,QR._1); // for (i <- 0 until n) curA(i)(i)+=Ann
          norm1 = norm2
          norm2 = normM(curA,null)
          iter-=1
        }
        if (iter == 0) throw new java.lang.IllegalArgumentException("Последовательность {Ak}, k=0:infinity не сходится по норме на "+ITER.toInt+" итерациях (с точностью до ERROR = "+ERROR+")")
        curA
      }
      // val xs = getSolve(QR_HR)
      def getSolve(QR: (Array[Array[Double]],Array[Array[Double]])): Array[Double] = {
        val Q = QR._1
        val R = QR._2
        val newB = new Array[Double](n)
        for (k <- 0 until n) newB(k) = MyVector.multi(b,Q,k)
        val xs = new Array[Double](n)
        // index
        var i = n-1
        while (i>=0) {
          Counter.countOpsMD+=2
          xs(i) = (1.0/R(i)(i))*(newB(i) - MyVector.multi(R(i),xs,i+1,n))
          i-=1;
        }
        xs
      }
      val eigenvalues = {          
        val diag = new Array[Double](n)
        for (i <- 0 until n) diag(i) = Ak(i)(i)
        diag
      }
    }
}
