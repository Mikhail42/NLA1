# NLA
Numerical Linear Algebra. 

GNU Affero General Public License

Original: 
В данной работе реализованы LU-,LUP-,QR- разложения, методы Якоби и Гаусса-Зейделя, методы нахождения решения СЛАУ через LU и LUP разложения, нахождения собственных чисел через QR-разложения. Дополнительно, реализованы многие операции с матрицами по подобию Matlab. 

К недостаткам работы относится долгое время выполнения: собственные числа на Scilab и Matlab вычисляются существенно быстрее. Для плохо обусловленной матрицы 200*200 время вычисления собственных чисел составляет порядка 30-40 секунд, против 2-3 секунд на Scilab (функция spec(A)). Также, не очень удачно реализован счетчик операций, -- впрочем, он хорошо выполняет свою роль. 


Translate: 
In this paper we implemented LU-, LUP-, QR- decomposition methods of Jacobi and Gauss-Seidel methods for finding solutions through Slough and LU LUP decomposition, finding eigenvalues ​​by QR-decomposition. Additionally, we implemented many matrix operations in the likeness of Matlab.

The disadvantages of the work is a long time: the eigenvalues ​​on Scilab and Matlab computed much faster. For ill-conditioned matrix 200 * 200 the computation time (eigenvalues) is about 30-40 seconds, 2-3 seconds against Scilab (spec(A)). Also, not very successfully implemented counter operations - however, it is well performed its role.

Bibliography:
 * [1] -- С. Ю. Гоголева, Матрицы и вычисления, 2004
 * [2] -- В. Б. Андреев, Численные методы
 * [3] -- J. E. Roman, The QR Algorithm. http://people.inf.ethz.ch/arbenz/ewp/Lnotes/chapter3.pdf
 * [4] -- G. Fasshauer  http://www.math.iit.edu/~fass/477577_Chapter_11.pdf
