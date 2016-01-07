package structLA

case class Matrix(A: Array[Array[Double]]) {
  val n = A.length
  def * (B: Array[Array[Double]]): Array[Array[Double]] = MyMatrix.multi(A, B)
  def mormMS = MyMatrix.normMS(A)
}