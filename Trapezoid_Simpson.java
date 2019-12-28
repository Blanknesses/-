import static java.lang.Math.*;

public class Trapezoid_Simpson {
    //区间数量
    private int N;

    public Trapezoid_Simpson(int n){
        N = n;
    }

    /**
     * 利用复化梯形求解e^3*cos(PI*x)在[a,b]上的积分
     * @param begin 起始值a
     * @param end 终止值b
     * @return 积分解
     */
    public double getTrapezoid(double begin, double end){
        double result = exp(3 * begin) * cos(PI * begin);
        double h = (end - begin) / N;
        for(int i = 1; i < N - 1; i ++){
            result = result + 2 * exp(3 * i * h) * cos(PI * i * h);
        }
        result = result + exp(3 * end) * cos(PI * end);
        result = result * h / 2;
        return result;
    }

    /**
     * 利用Simpson求解e^3*cos(PI*x)在[a,b]上的积分
     * @param begin 起始值a
     * @param end 终止值b
     * @return 积分解
     */
    public double getSimpson(double begin, double end){
        double result = exp(3 * begin) * cos(PI * begin) + exp(3 * end) * cos(PI * end);
        double h = (end - begin) / N;
        result = result + 4 * exp(3 * begin + h / 2) * cos(PI * (begin + h / 2));
        for(int i = 1; i < N - 1; i ++){
            result = result + 2 * exp(3 * i * h) * cos(PI * i * h) + 4 * exp(3 * i * h + h / 2) * cos(PI * (i * h + h / 2));
        }
        result = result * h / 6;
        return result;
    }
}
