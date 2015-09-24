package LR1

/**
 * @author misha
 * This object (!!!) is used in the calculation of chislel operations from the set {*, /, +, - =,?}.
 */
object Counter {
    /**  multi and delete */
    var countMD: Long = 0
    /**  plus, minus */
    var countPM: Long = 0
    /** Initialization numbers or assigning them values  */
    var countIA: Long = 0
    /** comparison  */
    var countC: Long = 0
    def setZeros(){countMD = 0; countPM = 0; countIA = 0; countC = 0;}
    override def toString() = (" {*,/} - "+countMD+"; {+,-} - "+countPM +"; {=} - "+countIA+"; {>,<,==} - "+countC)
}
