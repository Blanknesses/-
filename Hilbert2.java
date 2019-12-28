import java.text.DecimalFormat;

/**
 *  利用增广矩阵Gauss消元法实现
 *  误差向下积累
 */
public class Hilbert2 {
    private int N;
    private double[][] HilbertMatrix;
    private double[] result;
    private DecimalFormat df = new DecimalFormat( "0.000000 ");

    /**
     * 构造函数
     * 构造并初始化增广Hilbert矩阵，解向量
     * @param n
     */
    public Hilbert2(int n){
        N = n;
        HilbertMatrix = new double[n][n + 1];
        //Hilbert增广矩阵构造
        for(int i = 0; i < n; i++)
            for(int j = 0; j < n; j++)
                HilbertMatrix[i][j] = 1/(double)(i + j + 1);

        result = new double[n];
        for (int i = 0; i < n; i ++) {
            result[i] = 0;
            HilbertMatrix[i][n] = 1;
        }
    }

    /**
     * 增广矩阵按行进行Gauss消元
     */
    public void factorization(){
        int i = 0, j = 1;
        //按行进行消元
        for(; i < N - 1; i++, j++){
            for(int k = j; k < N; k ++){
                //获取行间商
                double Quotient = HilbertMatrix[k][i] / HilbertMatrix[i][i];
                if(Quotient == 0)
                    continue;
                //首位置0
                HilbertMatrix[k][i] = 0;
                //更新行
                if(Quotient != 0)
                    updateRaw(Quotient, k, j);
            }
        }
    }

    /**
     * 按行更新增广Hilbert矩阵
     * @param quotient 两行之间比例
     * @param row 当前待更新的行
     * @param column 待更新的行的第一个元素所在列
     */
    void updateRaw(double quotient, int row, int column){
        for(; column < N + 1; column ++)
            HilbertMatrix[row][column] = HilbertMatrix[row][column] - quotient * HilbertMatrix[row - 1][column];
    }

    /**
     * 直接逆序按行遍历增广矩阵，利用部分和化为1元1次方程组求解x
     */
    public void getResult(){
        double partSum = 0;
        //利用Hx = b求解x
        for(int i = N - 1; i >= 0; i --) {
            partSum = 0;
            //已知解计算的部分和，全部转化为一元一次方程求解
            for (int j = N - 1; j > i; j--) {
                partSum = partSum + result[j] * HilbertMatrix[i][j];
            }
            result[i] = (HilbertMatrix[i][N] - partSum) / HilbertMatrix[i][i];
        }
    }
    /**
     * 输出解向量结果
     */
    public void showResult(){
        for(int i = 0; i < N; i ++)
            System.out.println(df.format(result[i]));
    }

    /**
     * 输出矩阵结果
     */
    public void showMatrix(){
        for(int i = 0; i < N; i ++) {
            for (int j = 0; j < N + 1; j++)
                System.out.print(df.format(HilbertMatrix[i][j]) + " ");
            System.out.println();
        }
    }
}
