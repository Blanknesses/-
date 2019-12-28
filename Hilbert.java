import java.text.DecimalFormat;
import java.util.Arrays;

/**
 * 利用LU分解求解Hilbert矩阵
 * 误差向上积累
 */
public class Hilbert {
    //Hilbert矩阵的阶数
    private int N;
    private double[][] HilbertMatrix;
    private double[] b;
    private double[][] L_Matrix;
    private double[][] U_Matrix;
    private double[] result_x;
    private double[] result_y;
    private DecimalFormat df = new DecimalFormat( "0.000000 ");

    /**
     * 构造函数
     * 构造并初始化Hilbert矩阵，解向量，LU矩阵
     * @param n Hilbert矩阵的阶数
     */
    public Hilbert(int n){
        N = n;
        HilbertMatrix = new double[n][n];
        //Hilbert矩阵构造
        for(int i = 0; i < n; i++)
            for(int j = 0; j < n; j++)
                HilbertMatrix[i][j] = 1/(double)(i + j + 1);

        b = new double[n];
        result_x = new double[n];
        result_y = new double[n];
        for (int i = 0; i < n; i ++) {
            b[i] = 1;
            result_x[i] = 0;
            result_y[i] = 0;
        }

        L_Matrix = new double[n][n];
        U_Matrix = new double[n][n];
        //初始化LU矩阵
        for(int i = 0; i < n; i ++)
            for (int j = 0; j < n; j ++){
                L_Matrix[i][j] = 0;
                U_Matrix[i][j] = 0;
            }


    }

    /**
     * 对Hilbert矩阵进行LU分解
     * 更新LU矩阵，更新初始Hilbert矩阵将其化为上三角矩阵
     */
    public void LU_Factorization(){
        //初始化L矩阵对角元
        for(int i = 0; i < N; i ++)
            L_Matrix[i][i] = 1;
        int i = 0, j = 1;
        //按行进行消元
        for(; i < N - 1; i++, j++){
            //Hilbert矩阵不需要进行列主元消元
            //更新U矩阵第i行
            for(int k = i; k < N; k ++)
                U_Matrix[i][k] = HilbertMatrix[i][k];
            //更新L矩阵第i列
            for(int k = j; k < N; k ++){
                //获取行间商
                double Quotient = HilbertMatrix[k][i] / HilbertMatrix[i][i];
                if(Quotient == 0)
                    continue;
                //首位置0
                HilbertMatrix[k][i] = 0;
                L_Matrix[k][i] = Quotient;
                //更新行
                if(Quotient != 0)
                    updateRaw(Quotient, k, j);
            }
        }
        U_Matrix[N - 1][N - 1] = HilbertMatrix[N - 1][N - 1];
    }

    /**
     * LU分解，更新Hilbert矩阵的元素
     * @param quotient 两行之间比例
     * @param row 当前待更新的行
     * @param column 待更新的行的第一个元素所在列
     */
    void updateRaw(double quotient, int row, int column){
        for(; column < N; column ++)
            HilbertMatrix[row][column] = HilbertMatrix[row][column] - quotient * HilbertMatrix[row - 1][column];
    }

    /**
     * 计算解向量
     * 先根据Ly=b，计算y，正序遍历，利用部分和化为1元1次方程组求y
     * 再根据Ux=y，计算x，逆序遍历，同理
     */
    public void getResult(){
        //解向量初始化
        for(int i = 0; i < N; i ++) {
            result_y[i] = 0;
            result_x[i] = 0;
        }
        double partSum = 0;
        //利用Ly = b求解y
        for(int i = 0; i < N; i ++) {
            partSum = 0;
            //已知解计算的部分和，全部转化为一元一次方程求解
            for (int j = 0; j < i; j++) {
                partSum = partSum + result_y[j] * L_Matrix[i][j];
            }
            result_y[i] = (b[i] - partSum) / L_Matrix[i][i];
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

    /**
     * 输出矩阵结果
     */
    public void showMatrix() {
        for(int i = 0; i < N; i ++)
            System.out.println(Arrays.toString(HilbertMatrix[i]));
        System.out.println();
        for(int i = 0; i < N; i ++)
            System.out.println(Arrays.toString(L_Matrix[i]));
        System.out.println();
        for(int i = 0; i < N; i ++)
            System.out.println(Arrays.toString(U_Matrix[i]));
    }

    /**
     * 输出解向量结果
     */
    public void showResult(){
        System.out.println();
        for(int i = 0; i < N; i ++)
            System.out.println(df.format(result_y[i]));
        System.out.println();
        for(int i = 0; i < N; i ++)
            System.out.println(df.format(result_x[i]));
    }
}
