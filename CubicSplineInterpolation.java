import java.text.DecimalFormat;

import static java.lang.Math.pow;

public class CubicSplineInterpolation {
    //三次样条提供的点数
    private int N;
    //三次样条提供的点的值
    private double[] pointArr;
    //三次样条提供的点对应的值
    private double[] sArr;
    //严格对角占优系数矩阵
    private double[][] Matrix;
    //L阵
    private double[][] L_Matrix;
    //U阵
    private double[][] U_Matrix;
    //g向量
    private double[] g;
    //Ux = y解出的x解向量
    private double[] result_x;
    //Ly = b解出的y解向量
    private double[] result_y;

    private DecimalFormat df = new DecimalFormat( "0.000000 ");

    /**
     * 构造函数
     * 初始化各种矩阵
     * @param n 三次样条提供的点的值
     */
    CubicSplineInterpolation(int n){
        N = n;
        Matrix = new double[n][n];
        L_Matrix = new double[n][n];
        U_Matrix = new double[n][n];
        result_x = new double[n];
        result_y = new double[n];
        g = new double[n];
        for(int i = 0; i < n; i ++){
            g[i] = 0;
            result_y[i] = 0;
            result_x[i] = 0;
            for(int j = 0; j < n; j ++){
                if(i != j) {
                    Matrix[i][j] = 0;
                    L_Matrix[i][j] = 0;
                    U_Matrix[i][j] = 0;
                }
                else {
                    Matrix[i][j] = 2;
                    L_Matrix[i][j] = 1;
                }
            }
        }

        //已知题目为第二边界条件
        Matrix[0][1] = 1;
        Matrix[n - 1][n - 2] = 1;
    }

    /**
     * 计算lambda和mu，生成严格对角占优矩阵
     */
    public void getMatrix(){
        if(pointArr.length != 0){
            for(int i = 1; i < N - 1; i ++){
                Matrix[i][i - 1] = (pointArr[i + 1] - pointArr[i]) / ((pointArr[i + 1] - pointArr[i]) + (pointArr[i] - pointArr[i - 1]));
                Matrix[i][i + 1] = 1 - Matrix[i][i - 1];
            }
        }
    }
    /**
     * 对严格对角占优矩阵进行LU分解
     * 更新LU矩阵，更新初始严格对角占优矩阵将其化为上三角矩阵
     */
    public void LU_Factorization(){
        //初始化L矩阵对角元
        int i = 0, j = 1;
        //按行进行消元
        for(; i < N - 1; i++, j++){
            //Hilbert矩阵不需要进行列主元消元
            //更新U矩阵第i行
            for(int k = i; k < N; k ++)
                U_Matrix[i][k] = Matrix[i][k];
            //更新L矩阵第i列
            for(int k = j; k < N; k ++){
                //获取行间商
                double Quotient = Matrix[k][i] / Matrix[i][i];
                if(Quotient == 0)
                    continue;
                //首位置0
                Matrix[k][i] = 0;
                L_Matrix[k][i] = Quotient;
                //更新行
                if(Quotient != 0)
                    updateRaw(Quotient, k, j);
            }
        }
        U_Matrix[N - 1][N - 1] = Matrix[N - 1][N - 1];
    }

    /**
     * LU分解，更新严格对角占优矩阵的元素
     * @param quotient 两行之间比例
     * @param row 当前待更新的行
     * @param column 待更新的行的第一个元素所在列
     */
    void updateRaw(double quotient, int row, int column){
        for(; column < N; column ++)
            Matrix[row][column] = Matrix[row][column] - quotient * Matrix[row - 1][column];
    }

    /**
     * 根据边界条件值（此处为第二边界条件），计算常数向量
     * @param f0 第一个点处二阶导函数值
     * @param fn 第n个点处二阶导函数值
     */
    public void getBoundaryArr(double f0, double fn){
        g[0] = 3 * (sArr[1] - sArr[0]) / (pointArr[1] - pointArr[0]) - (pointArr[1] - pointArr[0]) / 2 * f0;
        g[N - 1] = 3 * (sArr[N - 1] - sArr[N - 2]) / (pointArr[N - 1] - pointArr[N - 2]) - (pointArr[N - 1] - pointArr[N - 2]) / 2 * f0;
        for(int i = 1; i < N - 1; i ++)
            g[i] = 3 * (Matrix[i][i + 1] * (sArr[i + 1] - sArr[i]) / (pointArr[i + 1] - pointArr[i]) +
                    Matrix[i][i - 1] * (sArr[i] - sArr[i - 1]) / (pointArr[i] - pointArr[i - 1]));
    }

    /**
     * 计算解向量
     * 先根据Ly=b，计算y，正序遍历，利用部分和化为1元1次方程组求y
     * 再根据Ux=y，计算x，逆序遍历，同理
     */
    public void getResult(){
        double partSum = 0;
        //利用Ly = b求解y
        for(int i = 0; i < N; i ++) {
            partSum = 0;
            //已知解计算的部分和，全部转化为一元一次方程求解
            for (int j = 0; j < i; j++) {
                partSum = partSum + result_y[j] * L_Matrix[i][j];
            }
            result_y[i] = (g[i] - partSum) / L_Matrix[i][i];
        }
        //利用Ux = y求解x
        for(int i = N - 1; i >= 0; i --) {
            partSum = 0;
            //已知解计算的部分和，全部转化为一元一次方程求解
            for (int j = N - 1; j > i; j--) {
                partSum = partSum + result_x[j] * U_Matrix[i][j];
            }
            result_x[i] = (result_y[i] - partSum) / U_Matrix[i][i];
        }
    }

    public void setPointArr(double[] pointArr) {
        this.pointArr = pointArr;
    }

    public void setsArr(double[] sArr){
        this.sArr = sArr;
    }

    /**
     * 输出对角占优矩阵
     */
    public void showMatrix(){
        for(int i = 0; i < N; i ++) {
            for (int j = 0; j < N; j++)
                System.out.print(df.format(Matrix[i][j]));
            System.out.println();
        }
        System.out.println();
        for(int i = 0; i < N; i ++) {
            for (int j = 0; j < N; j++)
                System.out.print(df.format(L_Matrix[i][j]));
            System.out.println();
        }
        System.out.println();
        for(int i = 0; i < N; i ++) {
            for (int j = 0; j < N; j++)
                System.out.print(df.format(U_Matrix[i][j]));
            System.out.println();
        }
    }

    public void showg(){
        for (int i = 0; i < N; i ++)
            System.out.println(g[i]);
    }


    public void showResult(){
        System.out.println();
        for(int i = 0; i < N; i ++)
            System.out.println(df.format(result_y[i]));
        System.out.println();
        for(int i = 0; i < N; i ++)
            System.out.println(df.format(result_x[i]));
    }

    public void getSx(){
        for(int i = 0; i < N - 1; i ++){
            String result = "区间[" + i + "," + (i + 1) + "]的方程为：";
            result = result + sArr[i] / (pointArr[i + 1] - pointArr[i]) * (pointArr[i + 1] - pointArr[i]) + " * (x - " + pointArr[i + 1] + ")^2 + "
                + df.format((2 * sArr[i] + (pointArr[i + 1] - pointArr[i]) * result_x[i]) / pow((pointArr[i + 1] - pointArr[i]),3)) + " * (x - " + pointArr[i] +")(x " + pointArr[i + 1] + ")^2 + "
                + sArr[i + 1] / ((pointArr[i + 1] - pointArr[i]) * (pointArr[i + 1] - pointArr[i])) + " * (x - " + pointArr[i] + ")^2 - "
                + df.format((2 * sArr[i + 1] + (pointArr[i + 1] - pointArr[i]) * result_x[i + 1]) / pow((pointArr[i + 1] - pointArr[i]),3)) + " * (x - " + pointArr[i + 1] +")(x " + pointArr[i] + ")^2";
            System.out.println(result);
        }
    }
}
