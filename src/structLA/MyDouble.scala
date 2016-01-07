package structLA;

import scala.math.abs
/**
 * @author Ionkin Mikhail
 */
object MyDouble {
  import math._
  
  var doubleError = 1e-10
  
  def compareDouble(a: Double, b: Double) = 
    if (abs(a-b)<doubleError) 0 else if (a>b) 1 else -1
  
  // !!! If you do otherwise, the program will work properly! 
  // Using in a QR-decomposition:  
  def sign(x: Double) = if (compareDouble(x, 0.0)>=0) 1 else -1
  
  def isZero(x: Double) = abs(x) < doubleError
}