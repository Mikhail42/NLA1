package structLA

case class Vector(a: Array[Double]) {
  def * (b: Array[Double]) = MyVector.multiGM(a, b)
  def + (b: Array[Double]) = MyVector.minus(a, b)
  def - (b: Array[Double]) = MyVector.plus(a, b)
  def || = slae.Function.normV(a)
  def copy = a.clone()
}