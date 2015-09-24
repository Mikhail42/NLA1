# NLA
Numerical Linear Algebra. 

Original: 
В данной работе реализованы LU-,LUP-,QR- разложения, методы Якоби и Гаусса-Зейделя, метод нахождения решения через LU и LUP разложения. 
Дополнительно, реализованы многие операции с матрицами по подобию Matlab. 
К недостаткам работы относится долгое время выполнения: собственные числа на Scilab и Matlab вычисляются существенно быстрее. Для плохо обусловленной матрицы 200*200 время вычислений составляет порядка 30-40 секунд, против 2-3 секунд на Scilab. Возможно, использование QR-алгоритма с двойными сдвигами и большая параллельность в вычислениях смогли бы существенно увеличить скорость. Также, не очень удачно реализован счетчик операций, -- впрочем, он хорошо выполняет свою роль. 

Translate: 
In this paper we implemented LU-, LUP-, QR- decomposition methods of Jacobi and Gauss-Seidel method for finding solutions through LUP and LU decomposition.
Additionally, we implemented many matrix operations in the likeness of Matlab.
The disadvantages of the work is a long time: the eigenvalues ​​on Scilab and Matlab computed much faster. For ill-conditioned matrix 200 * 200 the computation time is about 30-40 seconds, 2-3 seconds against Scilab. Perhaps the use of QR-algorithm with double shifts and more parallel to the calculations have been able to significantly increase the speed. Also, not very successfully implemented counter operations - however, it is well performed its role.

Bibliography:
 * [1] -- "С. Ю. Гоголева, Матрицы и вычисления, 2004
 * [2] -- "В. Б. Андреев, Численные методы"
 * [3] -- "J. E. Roman, The QR Algorithm. http://people.inf.ethz.ch/arbenz/ewp/Lnotes/chapter3.pdf"
 * [4] -- "G. Fasshauer  http://www.math.iit.edu/~fass/477577_Chapter_11.pdf"
