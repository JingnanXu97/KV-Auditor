package th_w;
import org.apache.commons.math3.distribution.BetaDistribution;

import java.util.Arrays;
import java.util.Random;
import java.util.ArrayList;
import java.util.List;
import org.apache.commons.math3.analysis.MultivariateFunction;
import org.apache.commons.math3.analysis.function.Exp;
import org.apache.commons.math3.analysis.function.Log;
import org.apache.commons.math3.analysis.function.Pow;
import org.apache.commons.math3.distribution.LaplaceDistribution;
import org.apache.commons.math3.optim.MaxEval;
import org.apache.commons.math3.optim.nonlinear.scalar.GoalType;
//import org.apache.commons.math3.optim.nonlinear.scalar.noderivative.BrentOptimizer;
//import org.apache.commons.math3.optim.nonlinear.scalar.noderivative.BrentOptimizer.Bracket;

public class the_ac {


    int iter = 100000000;
    int v1 = 0;
    int v2 = 1;
    int binNum = 4;
    double[] epsilon = {0.1, 0.2, 0.4, 0.8, 1.6, 3.2, 4, 6};


    public the_ac() throws Exception{
       audit();
    }


    public int getRandomChoice(double[] ue_val, int k, double thresh) {
        // Sum the elements of ue_val
        int sum = 0;
        for (double val : ue_val) {
            if(val > thresh)
                sum++;
        }

        Random random = new Random();

        // If the sum is 0, return a random integer between 0 and k-1
        if (sum == 0) {
            return random.nextInt(k);
        } else {
            // Create a list to store indices where ue_val[i] == 1
            List<Integer> indices = new ArrayList<>();
            for (int i = 0; i < ue_val.length; i++) {
                if (ue_val[i] > thresh) {
//                    System.out.println(ue_val[i]);
                    indices.add(i);
                }
            }


            int re =  indices.get(random.nextInt(indices.size()));
            return re;
        }
    }



    public void audit(){
        double[] thresh = {0.5124928545331653,0.5249429746156357,0.5495498751789415,0.596578396322232,0.6769874744974638,0.7823500118298539,0.8156771238097529 ,0.8691255209837462};
        //todo: the threshold is calculated with Python
        for(int eps = 0; eps < epsilon.length; eps++) {
            double[] return1 = new double[binNum];
            double[] return2 = new double[binNum];
            double emp_eps = 0.0;
            int count1 = 0;
            int count2 = 0;
            int test_statistic = 0;
            for (int i = 0; i < iter; i++) {
                return1 = THE(v1, -1, 1, epsilon[eps], thresh[eps]);
                int predict = getRandomChoice(return1, binNum, thresh[eps]);
                count1 += (predict == test_statistic) ? 1 : 0;
                return2 = THE(v2, -1, 1, epsilon[eps], thresh[eps]);

                predict = getRandomChoice(return2, binNum, thresh[eps]);
                count2 += (predict == test_statistic) ? 1 : 0;
            }

//            for(int i = 0; i < binNum; i++) {
                double[] cp1 = clopperPearson(count1, iter, 0.05);
                double[] cp2 = clopperPearson(count2, iter, 0.05);
                double log_m1 = Math.abs(Math.log(cp1[0] / cp2[1]));
                double log_m2 = Math.abs(Math.log(cp1[1] / cp2[0]));
                emp_eps = Math.abs(Math.max(log_m1, log_m2));
//            }
            System.out.println(emp_eps);



        }
    }

    //限制数据取值范围也是[-1, 1] bin的个数和kv中的设置相同
    public double[] cal_boundary(int b1, int b2, int binNum){
        double step = (double)(b2 - b1) / binNum;
        double[] boundary = new double[2 * binNum];

        boundary[0] = -1;
        boundary[binNum*2-1] = 1;

        for(int i = 1; i < binNum*2 - 1; i = i + 2){
            boundary[i] = boundary[i - 1] + step;
            boundary[i + 1] = boundary[i];
        }
        return boundary;
    }

    public static double laplace(double epsilon) {
        double sensitivity = 2.0;
        double scale = sensitivity / epsilon;
        LaplaceDistribution laplaceDist = new LaplaceDistribution(0, scale);
        return laplaceDist.sample();
    }


    // 对数据进行扰动
    public double[] perturbation(int value, double epsilon) {
        double[] encode = new double[binNum];
        // 编码，将特定位置置为 1
        for (int i = 0; i < binNum; i++) {
            encode[i] = (i == value) ? 1.0 : 0.0;
        }
        // 添加拉普拉斯噪声
        for (int i = 0; i < binNum; i++) {
            encode[i] = encode[i] + laplace(epsilon);
        }

        return encode;
    }

    public double[] THE(int value, int b1, int b2, double epsilon, double theta) {
        double[] perturb =  perturbation(value, epsilon);

        return perturb;  //support
    }

    public static double[] clopperPearson(int k, int n, double alpha) {
        double lowerBound, upperBound;


        BetaDistribution betaLower = new BetaDistribution(k, n - k + 1);
        lowerBound = betaLower.inverseCumulativeProbability(alpha / 2);


        BetaDistribution betaUpper = new BetaDistribution(k + 1, n - k);
        upperBound = betaUpper.inverseCumulativeProbability(1 - alpha / 2);
        return new double[]{lowerBound, upperBound};
    }
    public static void main(String[] args)  throws Exception {
        new the_ac();

    }

}

