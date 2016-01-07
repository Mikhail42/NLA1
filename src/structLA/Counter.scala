package structLA

/**
 * @author misha
 */
object Counter {
  
  /**  multi and delete */
  var countMD = 0L
  
  /**  plus, minus */
  var countPM = 0L
  
  /** Initialization numbers or assigning them values  */
  var countIA= 0L
  
  /** comparison  */
  var countC = 0L
  
  def setZeros = {countMD = 0;  countPM = 0; countIA = 0; countC = 0;}
  
  override def toString() = 
    (" {*,/} - "+countMD+"; {+,-} - "+countPM +"; {=} - "+countIA+"; {>,<,==} - "+countC)
}