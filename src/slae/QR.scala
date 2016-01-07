package slae

import structLA._
import math._

/**
 * @author Ionkin Mikhail
 * QR decomposition (also called a QR factorization) 

 * [1] -- "С. Ю. Гоголева, Матрицы и вычисления, 2004
 * [2] -- "В. Б. Андреев, Численные методы"
 * [3] -- "J. E. Roman, The QR Algorithm. http://people.inf.ethz.ch/arbenz/ewp/Lnotes/chapter3.pdf"
 * [4] -- "G. Fasshauer  http://www.math.iit.edu/~fass/477577_Chapter_11.pdf"
 */
class QR(A: Array[Array[Double]], b: Array[Double]) {  
  val n = A.length 
  lazy val QR = Function.QR.Decomposition.getQR_HR(A)
  lazy val eigenvalues = Function.QR.Eignvalues.eigenvalues(A)
  lazy val xs =  Function.QR.Solution.solve(QR, b)      
}