import static java.lang.Math.*;

public class Newton_Secant {
    private double[] coefficientArr;
    private int maxPower;
    private int iterationNum;

    /**
     * 构造函数
     * @param mp 多项式最高项系数
     */
    public Newton_Secant(int mp){
        maxPower = mp;
        coefficientArr = new double[mp + 1];
        iterationNum = 0;
    }

    /**
     * 输入系数向量，从低次到高次
     * @param coefficientArr 系数向量
     */
    public void setCoefficientArr(double[] coefficientArr) {
        this.coefficientArr = coefficientArr;
    }

    /**
     * 非线性方程组获取Newton迭代法的f(x)部分
     * @param x 当前x的值
     * @return f(x)的值
     */
    private double getFx(double x){
        double result = 0;
        for(int i = 0; i <= maxPower; i ++){
            if(coefficientArr[i] != 0)
                result = result + pow(x, i) * coefficientArr[i];
        }
        return result;
    }

    /**
     * 非线性方程组获取Newton迭代法的f'(x)部分
     * @param x 当前x的值
     * @return f‘(x)的值
     */
    private double getdF(double x){
        double result = 0;
        for(int i = 1; i <= maxPower; i ++){
            if(coefficientArr[i] != 0)
                result = result + i * pow(x, i - 1) * coefficientArr[i];
        }
        return result;
    }

    /**
     * 利用Newton法求非线性方程组
     * @param x 初始x值
     * @return Newton法迭代后获取的x值
     */
    public double Newton(double x){
        double x1 = x - getFx(x) / getdF(x);
        iterationNum = 1;
        while (abs(x1 - x) > 1.0E-6){
            x = x1;
            x1 = x - getFx(x) / getdF(x);
            iterationNum ++;
        }
        return x1;
    }

    /**
     * 非线性方程，自定义迭代次数Newton迭代法
     * @param x 初始x值
     * @param n 迭代次数
     * @return Newton法迭代后获取的x值
     */
    public double Newton(double x, int n){
        for(int i = 0; i < n; i ++)
            x = x - getFx(x) / getdF(x);
        return x;
    }

    /**
     * 利用割线法求非线性方程组
     * @param x0 初始x0的值
     * @param x1 初值x1的值
     * @return 割线法迭代后获取的x值
     */
    public double secantMethod(double x0, double x1){
        double x2 = x1 - getFx(x1) / (getFx(x1) - getFx(x0)) * (x1 - x0);
        iterationNum = 1;
        while (abs(x2 - x1) > 1.0E-6){
            x0 = x1;
            x1 = x2;
            x2 = x1 - getFx(x1) / (getFx(x1) - getFx(x0)) * (x1 - x0);
            iterationNum ++;
        }
        return x2;
    }

    /**
     * e^xsinx利用Newton法求解
     * @param x 初始x值
     * @return Newton法迭代后获取的x值
     */
    public double expSin_Newton(double x){
        double x1 = x - (exp(x) * sin(x)) / (exp(x) * (sin(x) + cos(x)));
        iterationNum = 1;
        while (abs(x1 - x) > 1.0E-6){
            x = x1;
            x1 = x - (exp(x) * sin(x)) / (exp(x) * (sin(x) + cos(x)));
            iterationNum ++;
        }
        return x1;
    }

    /**
     * e^xsinx利用割线法求解
     * @param x0 初始x0的值
     * @param x1 初值x1的值
     * @return 割线法迭代后获取的x值
     */
    public double expsin_SecantMethod(double x0, double x1){
        double x2 = x1 - (exp(x1) * sin(x1)) / (exp(x1) * sin(x1) - exp(x0) * sin(x0)) * (x1 - x0);
        iterationNum = 1;
        while (abs(x2 - x1) > 1.0E-6){
            x0 = x1;
            x1 = x2;
            x2 = x1 - (exp(x1) * sin(x1)) / (exp(x1) * sin(x1) - exp(x0) * sin(x0)) * (x1 - x0);
            iterationNum ++;
        }
        return x2;
    }

    public int getIterationNum(){
        return iterationNum;
    }
}
