package slae;

import slae._
import structLA._

/**
 * @author Ionkin Mikhail
 * [1] -- "С. Ю. Гоголева, Матрицы и вычисления, 2004
 * [2] -- "В. Б. Андреев, Численные методы"
 * [3] -- "J. E. Roman, The QR Algorithm. http://people.inf.ethz.ch/arbenz/ewp/Lnotes/chapter3.pdf"
 * [4] -- "G. Fasshauer  http://www.math.iit.edu/~fass/477577_Chapter_11.pdf"
 */
case class SLAE(A : Array[Array[Double]], x: Array[Double]){
  val n = A.length
  lazy val b = MyMatrix.multi(A,x)
  lazy val lu  = new LU(A,b)
  lazy val lup = new LUP(A,b)
  lazy val qr  = new QR(A,b)
  lazy val jacobiAndSeidel = new JacobiAndSeidel(A,b) 
}