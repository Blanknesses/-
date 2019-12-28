import java.text.DecimalFormat;

import static java.lang.Math.abs;

public class Jacobi_GS {
    private int N;
    //记录迭代次数
    private int iterationNum;
    private double[][] Matrix;
    private double[] b;
    private double[] result;
    private DecimalFormat df = new DecimalFormat( "0.000000 ");

    /**
     * 构造并初始化严格对角占优矩阵，解向量，常数向量
     * @param n 严格对角占优矩阵阶数
     */
    public Jacobi_GS(int n){
        N = n;
        Matrix = new double[n][n];
        b = new double[n];
        result = new double[n];
        //构造并初始化严格对角占优矩阵，解向量，常数向量
        for(int i = 0; i < n; i ++) {
            for (int j = 0; j < n; j++) {
                Matrix[i][j] = 0;
            }
            result[i] = 0;
            b[i] = 1;
        }

        //严格对角占优矩阵，常数向量赋值
        b[0] = 2;
        b[n - 1] = 2;
        Matrix[0][0] = 3;
        Matrix[0][1] = -1;
        Matrix[n - 1][n - 2] = -1;
        Matrix[n - 1][n - 1] = 3;
        for(int i = 1; i < n - 1; i ++){
                Matrix[i][i] = 3;
                Matrix[i][i - 1] = -1;
                Matrix[i][i + 1] = -1;
        }

        iterationNum = 0;
    }

    /**
     * 进行Jacobi n次迭代
     */
    public void Jacobi(){
        while (getJacobi())
            iterationNum ++;
    }

    /**
     * 进行Jacobi迭代获取解向量
     */
    private boolean getJacobi(){
        //因为Jacobi迭代完成一次迭代后才更新x，利用temp数组暂存更新的x
        double[] temp = new double[N];
        for(int i = 0; i < N; i ++){
            //计算系数矩阵当前行与解向量乘积的部分和
            double partSum = 0;
            for(int j = 0; j < N; j ++)
                partSum = partSum + Matrix[i][j] * result[j];
            //将对角元素排除
            partSum = partSum - Matrix[i][i] * result[i];
            temp[i] = 1 / Matrix[i][i] * (b[i] - partSum);
        }
        //判断解向量精度
        boolean ret = judgeResult(temp);
        //更新解向量
        for(int i = 0; i < N; i ++)
            result[i] = temp[i];
        return ret;
    }

    /**
     * 进行GS n次迭代
     */
    public void GS(){
        while (getGS())
            iterationNum ++;
    }

    /**
     * 进行GS迭代获取解向量
     */
    private Boolean getGS(){
        double[] temp = new double[N];
        for(int i = 0; i < N; i ++)
            temp[i] = result[i];
        for(int i = 0; i < N; i ++){
            //计算系数矩阵当前行与解向量乘积的部分和
            double partSum = 0;
            for(int j = 0; j < N; j ++)
                partSum = partSum + Matrix[i][j] * result[j];
            //将对角元素排除
            partSum = partSum - Matrix[i][i] * result[i];
            result[i] = 1 / Matrix[i][i] * (b[i] - partSum);
        }
        //判断解向量精度
        boolean ret = judgeResult(temp);
        return ret;
    }

    /**
     * 比较两次迭代解之间的变化量，当解向量的变化量小于1.0E-6时停止迭代
     * @param temp Jacobi迭代的第i次迭代解向量，GS迭代的第i-1次迭代向量
     * @return 不满足精度，返回true，迭代继续，满足返回false，结束迭代
     */
    private boolean judgeResult(double[] temp){
        for(int i = 0; i < N; i ++)
            if(abs(temp[i] - result[i]) > 1.0E-6)
                return true;
        return false;
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
    }

    /**
     * 输出解向量
     */
    public void showResult(){
        for(int i = 0; i < N; i ++)
            System.out.print(df.format(result[i]));
        System.out.println();
    }

    public int getIterationNum(){
        return iterationNum;
    }


}
