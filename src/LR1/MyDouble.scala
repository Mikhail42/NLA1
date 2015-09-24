package LR1;
/**
 * @author Ionkin Mikhail
 */
object MyDouble {
  import math._
  val ERROR = 1e-8
  def compareDouble(a: Double, b: Double) = if (abs(a-b)<ERROR) 0 else if (a>b) 1 else -1
  // !!! If you do otherwise, the program will work properly!
  def sign(x: Double) = if (compareDouble(x,0.0)>=0) 1 else -1
  def isZero(x: Double) = abs(x)<ERROR
}
