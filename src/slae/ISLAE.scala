package slae
/**
 * @author Ionkin Mikhail
 * some static methods for SLAE

 * [1] -- "С. Ю. Гоголева, Матрицы и вычисления, 2004
 * [2] -- "В. Б. Андреев, Численные методы"
 * [3] -- "J. E. Roman, The QR Algorithm. http://people.inf.ethz.ch/arbenz/ewp/Lnotes/chapter3.pdf"
 * [4] -- "G. Fasshauer  http://www.math.iit.edu/~fass/477577_Chapter_11.pdf"
 */
object ISLAE {
    import LR1._
     /** Jacobi or Seidel 
     *  1 - Jacobi method
     *  2 - Seidel method
     *  */
    var METHOD_ID = 1
    /** 
     *  p of p-norms: p = NORM_ID
     *  if (p==3) => (p = Infinity)
     *  */
    var NORM_ID = 3
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
    /**
     * LU-decomposition. 
     * [1 : 1.1,1.2], [2: 1.2, 4.4] 
     * @param dwpFlag -- (true => decomposition with pivoting, false => simple decomposition)
     * @q -- the resulting vector permutations (if dwpFlag == true). 
     * @return (L,U). q is UPDATE
     */
    @throws(classOf[java.lang.IllegalArgumentException])
    def getLU(matr: Array[Array[Double]], dwpFlag: Boolean, q: Array[Int]): (Array[Array[Double]],Array[Array[Double]]) = {       
      val n = matr.length
      val U = Array.ofDim[Double](n,n)
      val L = Array.ofDim[Double](n,n)
      Counter.countIA+=n*n; Counter.countC+=n*n+n; Counter.countPM+=n*n
      /**
       *  [2]: 
       *   стр. 11: Можно чередовать вычисление строк U(k)(:) и столбцов L(:)(k)
       *   стр. 42: Устойчивость не гарантирована 
       */
      for (k <- 0 until n){
        for (j <- k until n) U(k)(j) = matr(k)(j) - MyVector.multi(L(k),U,j,0,k);
        if (dwpFlag){
          // индекс элемента с наибольшим модулем
          val ind = MyVector.indexOfMaxAbs(U(k),k)
          MyMatrix.swapColumn(matr,k,ind); MyMatrix.swapColumn(U,k,ind); 
          Counter.countC+=1;
          if (k!=ind) {Counter.countIA+=3; val g = q(k); q(k)=q(ind); q(ind) = g;}
        }
        if (MyDouble.isZero(U(k)(k))) throw new java.lang.IllegalArgumentException("LU-разложение невозможно ввиду необходимости деления на сверхмалое число")
        Counter.countMD+=(n-k); Counter.countIA+=(n-k+1); Counter.countPM+=(n-k); Counter.countC+=(n-k)
        val m1Ukk = 1.0/U(k)(k); L(k)(k) = 1.0;
        for (i <- k+1 until n) L(i)(k) = (matr(i)(k) - MyVector.multi(L(i),U,k,0,k))*m1Ukk; 
      }
      (L,U)
    }
    /**
    * [1], p. 7, eq. (1.5)
    */
    def solve(L: Array[Array[Double]], U: Array[Array[Double]], b: Array[Double]): Array[Double] = {
      val n = L.length
      Counter.countMD+=n; Counter.countIA+=(n<<2); Counter.countPM+=(n+n); Counter.countC+=(n+n)
      val ys = new Array[Double](n)
      import LR1.MyVector.multi
      for (i <- 0 until n) ys(i) = b(i) - multi(L(i),ys,0,i)
      val xs = new Array[Double](n)
      var i = n-1
      while (i>=0) {
        xs(i) = (1.0/U(i)(i))*(ys(i) - multi(U(i),xs,i+1,n))
        i-=1
      }
      xs
    }
}
