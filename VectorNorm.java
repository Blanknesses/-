import static java.lang.Math.abs;
import static java.lang.Math.sqrt;

public class VectorNorm {
    private double[] vector;
    private int N;

    /**
     * 向量构造函数1
     * 构造成1,1/2,...,1/n的向量
     * @param n 向量长度
     */
    private void generateVector1(int n){
        vector = new double[n];
        for (int i = 1; i <= n; i++)
            vector[i - 1] = 1/(double)i;
    }

    /**
     * 向量构造函数2
     * 构造成1,2,...,n的向量
     * @param n 向量长度
     */
    private void generateVector2(int n){
        vector = new double[n];
        for (int i = 1; i <= n; i++)
            vector[i - 1] = i;
    }

    /**
     * 构造函数，调用向量构造函数
     * @param n 向量长度
     */
    public VectorNorm(int n){
        N = n;
        //generateVector1(n);
        generateVector2(n);
    }

    /**
     * 向量1-范数计算
     * @return 向量元素和
     */
    public double norm1(){
        double result = 0;
        for(int i = 0; i < N; i ++)
            result = result + vector[i];
        return result;
    }

    /**
     * 向量2-范数计算
     * @return 向量元素平方和开方
     */
    public double norm2(){
        double result = 0;
        for(int i = 0; i < N; i ++)
            result = result + vector[i] * vector[i];
        return sqrt(result);
    }

    /**
     * 向量无穷范数计算
     * @return 向量最大元素
     */
    public double normInfty(){
        double max = vector[0];
        for(int i = 1; i < N; i ++)
            if(vector[i] > max)
                max = vector[i];
        return max;
    }

    public double[] getVector() {
        return vector;
    }
}
