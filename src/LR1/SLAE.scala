package LR1;
/**
 * @author Ionkin Mikhail
 * [1] -- "С. Ю. Гоголева, Матрицы и вычисления, 2004
 * [2] -- "В. Б. Андреев, Численные методы"
 * [3] -- "J. E. Roman, The QR Algorithm. http://people.inf.ethz.ch/arbenz/ewp/Lnotes/chapter3.pdf"
 * [4] -- "G. Fasshauer  http://www.math.iit.edu/~fass/477577_Chapter_11.pdf"
 */
class SLAE(A : Array[Array[Double]], x: Array[Double]){
  import slae._
  val b = MyMatrix.multi(A,x)
  
  private var lu: LU = null;  def getLU = lu; def setLU={lu = new LU(A,b)}
  private var lup: LUP = null; def getLUP = lup; def setLUP={lup = new LUP(A,b)}
  private var qr: QR = null; def getQR = qr; def setQR={qr =  new slae.QR(A,b)}
  private var jacobiAndSeidel: JacobiAndSeidel = null; def getJacobiAndSeidel = jacobiAndSeidel; def setJacobiAndSeidel={jacobiAndSeidel = new JacobiAndSeidel(A,b)}
}
