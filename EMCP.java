import java.util.Arrays;

import static java.lang.Math.abs;

public class EMCP {
    //定义行列
    private static int Row = 9;
    private static int Column = 9;
    //定义两个解向量
    private double[] result_x = new double[Column];
    private double[] result_y = new double[Column];
    //列主元LU分解，置换矩阵P
    private double[][] P;
    //系数矩阵
    private double[][] Matrix = {
            {31,-13,0,0,0,-10,0,0,0},
            {-13,35,-9,0,-11,0,0,0,0},
            {0,-9,31,-10,0,0,0,0,0},
            {0,0,-10,79,-30,0,0,0,-9},
            {0,0,0,-30,57,-7,0,-5,0},
            {0,0,0,0,-7,47,-30,0,0},
            {0,0,0,0,0,-30,41,0,0},
            {0,0,0,0,-5,0,0,27,-2},
            {0,0,0,-9,0,0,0,-2,29}};
    //解向量
    private double[] b = {-15,27,-23,0,-20,12,-7,7,10};

    //LU矩阵
    private double[][] L_Matrix = new double[Row][Column];
    private double[][] U_Matrix = new double[Row][Column];

    //逆矩阵
    private double[][] inverseMatrix = new double[Row][Column];

    public EMCP(){
        P = new double[Row][Column];
        for(int i = 0; i < Row; i ++)
            for(int j = 0; j < Column; j ++){
                if(i == j)
                    P[i][j] = 1;
                else
                    P[i][j] = 0;
                inverseMatrix[i][j] = 0;
            }
    }

    /**
     * 对矩阵进行列主元LU分解
     * 更新LU矩阵，更新初始矩阵将其化为上三角矩阵
     */
    public void LU_cpFactorization(){
        //初始化L矩阵对角元
        for(int i = 0; i < Row; i ++)
            L_Matrix[i][i] = 1;
        int i = 0, j = 1;
        double Max;
        //按行进行消元
        for(; i < Row - 1; i++, j++){
            Max = Matrix[i][i];
            int target = i;
            //找出最大列主元所在行
            for(int k = j; k < Column; k++){
                if(abs(Max) < abs(Matrix[k][i])) {
                    Max = Matrix[k][i];
                    target = k;
                }
            }
            //若初始行不是最大列主元所在行，进行行变换
            if(i != target)
                swap(i, target);
            //更新U矩阵第i行
            for(int k = i; k < Column; k ++)
                U_Matrix[i][k] = Matrix[i][k];
            //更新L矩阵第i列
            for(int k = j; k < Row; k ++){
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
        U_Matrix[Row - 1][Column - 1] = Matrix[Row - 1][Column - 1];
    }

    /**
     * 对矩阵进行LU分解
     * 更新LU矩阵，更新初始矩阵将其化为上三角矩阵
     */
    public void LU_Factorization(){
        //初始化L矩阵对角元
        for(int i = 0; i < Row; i ++)
            L_Matrix[i][i] = 1;
        int i = 0, j = 1;
        double Max;
        //按行进行消元
        for(; i < Row - 1; i++, j++){
            //更新U矩阵第i行
            for(int k = i; k < Column; k ++)
                U_Matrix[i][k] = Matrix[i][k];
            //更新L矩阵第i列
            for(int k = j; k < Row; k ++){
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
        U_Matrix[Row - 1][Column - 1] = Matrix[Row - 1][Column - 1];
    }

    /**
     * 寻找列主元后进行行变换
     * @param i 待交换行
     * @param j 目标交换行
     */
    public void swap(int i, int j){
        double temp = 0;
        for (int k = 0; k < Column; k++){
            temp = Matrix[i][k];
            Matrix[i][k] = Matrix[j][k];
            Matrix[j][k] = temp;

            temp = P[i][k];
            P[i][k] = P[j][k];
            P[j][k] = temp;
        }
    }

    /**
     * LU分解，更新矩阵的元素
     * @param quotient 两行之间比例
     * @param row 当前待更新的行
     * @param column 待更新的行的第一个元素所在列
     */
    void updateRaw(double quotient, int row, int column){
        for(; column < Column; column ++)
            Matrix[row][column] = Matrix[row][column] - quotient * Matrix[row - 1][column];
    }

    /**
     * 计算解向量
     * 先根据Ly=b，计算y，正序遍历，利用部分和化为1元1次方程组求y
     * 再根据Ux=y，计算x，逆序遍历，同理
     */
    public void getResult(){
        //解向量初始化
        for(int i = 0; i < Column; i ++) {
            result_y[i] = 0;
            result_x[i] = 0;
        }
        double partSum = 0;
        //利用Ly = b求解y
        for(int i = 0; i < Column; i ++) {
            partSum = 0;
            //已知解计算的部分和，全部转化为一元一次方程求解
            for (int j = 0; j < i; j++) {
                partSum = partSum + result_y[j] * L_Matrix[i][j];
            }
            result_y[i] = (b[i] - partSum) / L_Matrix[i][i];
        }
        //利用Ux = y求解x
        for(int i = Column - 1; i >= 0; i --) {
            partSum = 0;
            //已知解计算的部分和，全部转化为一元一次方程求解
            for (int j = Column - 1; j > i; j--) {
                partSum = partSum + result_x[j] * U_Matrix[i][j];
            }
            result_x[i] = (result_y[i] - partSum) / U_Matrix[i][i];
        }
    }

    /**
     * 计算原矩阵的逆
     * 利用Ly = ei Uai = y
     */
    public void getInverseMatrix(){
        double partSum = 0;
        for(int k = 0; k < Row; k ++) {
            //利用Ly = b求解y
            for(int i = 0; i < Row; i ++)
                b[i] = 0;
            for(int i = 0; i < Column; i ++) {
                result_y[i] = 0;
                result_x[i] = 0;
            }
            b[k] = 1;
            for (int i = 0; i < Column; i++) {
                partSum = 0;
                //已知解计算的部分和，全部转化为一元一次方程求解
                for (int j = 0; j < i; j++) {
                    partSum = partSum + result_y[j] * L_Matrix[i][j];
                }
                result_y[i] = (b[i] - partSum) / L_Matrix[i][i];
            }
            //利用Ux = y求解x
            for (int i = Column - 1; i >= 0; i--) {
                partSum = 0;
                //已知解计算的部分和，全部转化为一元一次方程求解
                for (int j = Column - 1; j > i; j--) {
                    partSum = partSum + result_x[j] * U_Matrix[i][j];
                }
                result_x[i] = (result_y[i] - partSum) / U_Matrix[i][i];
            }

            for(int i = 0; i < Row; i ++){
                inverseMatrix[i][k] = result_x[i];
            }
        }
    }

    /**
     * 因为LU矩阵分别为下，上三角阵，行列式为对角元之积
     * @return 行列式的值
     */
    public double getDeterminant(){
        double l = 1;
        double u = 1;
        for(int i = 0; i < Row; i ++){
            l = l * L_Matrix[i][i];
            u = u * U_Matrix[i][i];
        }
        return l * u;
    }

    public void showMatrix() {
        for(int i = 0; i < 9; i ++)
            System.out.println(Arrays.toString(Matrix[i]));
        System.out.println();
        for(int i = 0; i < 9; i ++)
            System.out.println(Arrays.toString(L_Matrix[i]));
        System.out.println();
        for(int i = 0; i < 9; i ++)
            System.out.println(Arrays.toString(U_Matrix[i]));
        System.out.println();
        for(int i = 0; i < 9; i ++)
            System.out.println(Arrays.toString(inverseMatrix[i]));
    }

    public void showResult(){
        System.out.println();
        System.out.println(Arrays.toString(result_y));
        System.out.println();
        System.out.println(Arrays.toString(result_x));
    }

    public static void main(String[] args) {
        EMCP test = new EMCP();
        test.LU_cpFactorization();
        test.getResult();
        test.showResult();
        System.out.println(test.getDeterminant());
    }
}
